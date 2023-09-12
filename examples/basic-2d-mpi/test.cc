#include <mpi.h>
#include <iostream>
#include <memory>
#include <math.h>

#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/mpi/gallery.h>
#include <cedar/2d/mpi/tausch_exchanger.h>

#include <cedar/util/time_log.h>
#include "vcycle.hpp"

#include <ftl/KernelRegistry.hpp>
#include <ftl/Buffer.hpp>
#include <ftl/Runtime.hpp>
#include <ftl/String.hpp>
#include <ftl/Device.hpp>


static void set_problem(cedar::cdr2::mpi::grid_func & b) {
    using namespace cedar;
    using namespace cedar::cdr2;

    const double pi = M_PI;

    auto rhs = [pi](real_t x, real_t y) {
        return 8*(pi*pi)*sin(2*pi*x)*sin(2*pi*y);
    };

    auto & topo = b.grid();

    b.set(0);

    real_t igs = topo.is(0);
    real_t jgs = topo.is(1);

    real_t hx = 1.0 / (topo.nglobal(0) - 1);
    real_t hy = 1.0 / (topo.nglobal(1) - 1);

    real_t h2 = hx*hy;

    real_t nlx = topo.nlocal(0) - 2;
    real_t nly = topo.nlocal(1) - 2;

    real_t i1 = nlx + 1;
    real_t j1 = nly + 1;

    for (auto j : range<len_t>(1, j1)) {
        for (auto i : range<len_t>(1, i1)) {
            len_t is = igs + i;
            len_t js = jgs + j;

            real_t x = (is-1)*hx;
            real_t y = (js-1)*hy;

            b(i,j) = rhs(x, y) * h2;

        }
    }
}


static void set_solution(cedar::cdr2::mpi::grid_func & q) {
    using namespace cedar;

    const double pi = M_PI;

    auto sol = [pi](real_t x, real_t y) {
        return sin(2*pi*x)*sin(2*pi*y);
    };

    auto & topo = q.grid();

    real_t igs = topo.is(0);
    real_t jgs = topo.is(1);

    real_t hx = 1.0 / (topo.nglobal(0) - 1);
    real_t hy = 1.0 / (topo.nglobal(1) - 1);

    for (auto j : q.range(1)) {
        for (auto i : q.range(0)) {
            len_t is = igs + i;
            len_t js = jgs + j;

            real_t x = (is-1)*hx;
            real_t y = (js-1)*hy;

            q(i,j) = sol(x,y);
        }
    }
}


template <typename dtype, unsigned short ND>
ftl::Buffer<dtype> to_buffer(cedar::array<dtype, ND>& array) {
    std::vector<uint32_t> shape(ND);

    for (std::size_t i = 0; i < ND; ++i) {
        const uint32_t dim = array.len(i);
        shape[i] = dim;
    }

    ftl::Buffer<dtype> buf(array.data(), shape);

    return buf;
}


template <typename len_t, typename dtype, unsigned short ND>
ftl::Buffer<dtype> to_buffer(cedar::aarray<len_t, dtype, ND>& array) {
    std::vector<uint32_t> shape(ND);

    for (std::size_t i = 0; i < ND; ++i) {
        const uint32_t dim = array.len(i);
        shape[i] = dim;
    }

    ftl::Buffer<dtype> buf(array.data(), shape);

    return buf;
}


template <typename stencil_type>
ftl::Buffer<real_t> so_to_buffer(cedar::cdr2::mpi::stencil_op<stencil_type>& array) {
    std::vector<uint32_t> shape = {array.len(0) + 1, array.len(1) + 1, array.len(2)};
    ftl::Buffer<real_t> buf(array.data(), shape);

    return buf;
}


template <typename T>
void save_file(ftl::Buffer<T>& buf, const std::string& fname) {
    std::ofstream file(fname);
    file << buf.tostr_unformatted() << std::endl;
    file.close();
}


void test_halo(long nx, long ny,
               int rank, int size,
               std::vector<cedar::topo_ptr>& topos,
               ftl::Buffer<real_t>& so,
               ftl::Buffer<real_t>& b,
               cedar::cdr2::mpi::halo_exchange* halof) {

    ftl::Buffer<real_t> x({nx, ny});
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            x(i, j) = static_cast<real_t>(rank + 1);
        }
    }

    std::cout << "Before halo exchange: " << std::endl;
    std::cout << x << std::endl;

    halof->exchange_func(1, x.data());

    std::cout << "After halo exchange: " << std::endl;
    std::cout << x << std::endl;
}


void test_setup_interp(
    long nx, long ny,
    int rank, int size,
    std::vector<cedar::topo_ptr>& topos,
    ftl::Buffer<real_t>& so,
    ftl::Buffer<real_t>& b,
    cedar::cdr2::mpi::halo_exchange* halof) {

    ftl::Buffer<real_t> soc({topos[1]->nlocal(0) + 1, topos[1]->nlocal(1) + 1, 5}, ftl::Buffer_Both);
    ftl::Buffer<real_t> ci({topos[1]->nlocal(0), topos[1]->nlocal(1), 8}, ftl::Buffer_Both);

    int kc = 1;
    int kf = 2;
    int nog = 2;

    cedar::grid_topo::igrd_t igrd = topos[0]->get_igrd();
    ftl::Buffer<uint32_t> igrd_buf(igrd->data(), {nog, NBMG_pIGRD});

    int nstncl = 3;
    int jpn = BMG_BCs_definite;
    int ifd = 1;

    void* halof_void = static_cast<void*>(halof);

    std::cout <<  "Coarsening from " <<
        topos[0]->nlocal(0) << " x " << topos[0]->nlocal(1) << " to " <<
        topos[1]->nlocal(0) << " x " << topos[1]->nlocal(1) << std::endl;

    MPI_BMG2_SymStd_SETUP_interp_OI<ftl::device::GPU>(
        kf, kc, so, ci,
        topos[0]->nlocal(0), topos[0]->nlocal(1),
        topos[1]->nlocal(0), topos[1]->nlocal(1),
        nog, nog, igrd_buf, ifd, nstncl, jpn, halof_void);

    int local_i = topos[0]->nlocal(0);
    int local_j = topos[0]->nlocal(1);
    int local_ic = topos[1]->nlocal(0);
    int local_jc = topos[1]->nlocal(1);

    std::cout << "MPI_BMG2_SymStd_SETUP_interp_OI(" << kf << ", " << kc << ", fopd, P, " <<
        local_i << ", " << local_j << ", " <<
        local_ic << ", " << local_jc << ", " <<
        nog << ", " << nog << ", igrd, " << ifd << ", " << nstncl << ", " << jpn <<
        ", halof)" << std::endl;

    std::cout << std::endl;
    std::cout << ci << std::endl;
}


int main(int argc, char *argv[]) {
    ftl::runtime_parse_args(argc, argv);
    ftl::device::autodetect();
    ftl::load_kernels(true);

    using namespace cedar;
    using namespace cedar::cdr2;

    int provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    /* Put a barrier here in case we are debugging a single process */
    MPI_Barrier(MPI_COMM_WORLD);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cout << "[Rank " << rank << "]: " << "Running on " << ftl::device::get_name() << std::endl;

    timer_init(MPI_COMM_WORLD);

    config conf;
    auto ndofs = conf.getvec<len_t>("grid.n");
    auto fgrid = util::create_topo(conf);
    fgrid->comm = MPI_COMM_WORLD;
    fgrid->grow(2);

    auto cgrid = std::make_shared<cedar::grid_topo>(fgrid->get_igrd(), 0, 2);

    /* Create the coarse grid topology.  This is stolen from Cedar */
    {
        len_t NLxg = fgrid->nlocal(0) - 2;
        len_t NLyg = fgrid->nlocal(1) - 2;
        len_t NGxg = (fgrid->nglobal(0) - 1) / 2 + 2;
        len_t NGyg = (fgrid->nglobal(1) - 1) / 2 + 2;

        cgrid->nglobal(0) = NGxg;
        cgrid->nglobal(1) = NGyg;

        if ((fgrid->is(0) % 2) == 1) {
            cgrid->is(0) = (fgrid->is(0) + 1) / 2;
            NLxg = (NLxg + 1) / 2;
        } else {
            cgrid->is(0) = fgrid->is(0)/2 + 1;
            if (NLxg % 2 == 1) {
                NLxg = (NLxg-1)/2;
            } else {
                NLxg = (NLxg+1)/2;
            }
        }

        if (fgrid->is(1) % 2 == 1) {
            cgrid->is(1) = (fgrid->is(1)+1) / 2;
            NLyg = (NLyg+1) / 2;
        } else {
            cgrid->is(1) = fgrid->is(1) / 2 + 1;
            if (NLyg % 2 == 1) {
                NLyg = (NLyg - 1) / 2;
            } else {
                NLyg = (NLyg+1)/2;
            }
        }

        cgrid->nlocal(0) = NLxg + 2;
        cgrid->nlocal(1) = NLyg + 2;

        cgrid->nproc(0) = fgrid->nproc(0);
        cgrid->nproc(1) = fgrid->nproc(1);
        cgrid->coord(0) = fgrid->coord(0);
        cgrid->coord(1) = fgrid->coord(1);
    }

    /* Set up halo exhange using the legacy MSG exchanger */
    /* This should set up the intermediate data structures when setup is called...? */
    auto kman = cedar::cdr2::mpi::build_kernel_manager(conf);
    auto sman = kman->services();
    auto halo_exchange = sman.fortran_handle<cedar::cdr2::mpi::halo_exchange>();
    auto& halo_service = sman.get<cedar::cdr2::mpi::halo_exchange>();
    auto tausch_service = dynamic_cast<cedar::cdr2::mpi::tausch_exchanger*>(&halo_service);
    if (tausch_service == nullptr) {
        throw std::runtime_error("tausch service is null");
    }
    std::vector<cedar::topo_ptr> topos;
    topos.push_back(fgrid);
    topos.push_back(cgrid);
    halo_service.setup(topos);

    auto nx = fgrid->nlocal(0);
    auto ny = fgrid->nlocal(1);

    /* Create the stencil operator */
    auto so_cedar = mpi::gallery::poisson(fgrid);
    ftl::Buffer<real_t> so = so_to_buffer<five_pt>(so_cedar);
    int kf = 2;
    // auto ctx = tausch_service->context();

    // std::cout << "ctx msg_geom (" << ctx.msg_geom.size() << ")" << std::endl;
    // auto msg_geom_ptr = ctx.msg_geom.data();
    // for (std::size_t i = 0; i < ctx.msg_geom.size(); ++i) {
    //     std::cout << msg_geom_ptr[i] << " ";
    // }
    // std::cout << std::endl;

    /* On rank zero, we get values like...
       0 0 0 0 0 0 0 0 0 0 4294967271 26 27 28 29 30 4294967265 32 33 34 35 36 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                           ^ -27                     ^ -31

                           Are those large/negative values supposed to be in there...?
                           On rank one this is causing indexing out of bounds issues.
    */

    // I assume we need to do a halo exchange on the fine grid stencil?
    // int fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
    // BMG2_SymStd_SETUP_fine_stencil(kf, so, nx, ny, 3,
    //                                ftl::Buffer<unsigned int>(ctx.msg_geom.data(), ctx.msg_geom.size()), ctx.msg_geom.size(),
    //                                to_buffer(ctx.pMSGSO),
    //                                ftl::Buffer<real_t>(ctx.msg_buffer.data(), ctx.msg_buffer.size()), ctx.msg_buffer.size(),
    //                                fcomm);


    tausch_service->run(so_cedar);

    /* Setup RHS */
    mpi::grid_func b_cedar(fgrid);
    set_problem(b_cedar);
    ftl::Buffer<real_t> b = to_buffer(b_cedar);

    /* Setup exact solution */
    mpi::grid_func exact_sol_cedar(fgrid);
    set_solution(exact_sol_cedar);
    ftl::Buffer<real_t> exact_sol = to_buffer(exact_sol_cedar);



    //test_halo(nx, ny, rank, size, fgrid, so, b, halo_exchange);
    test_setup_interp(nx, ny, rank, size, topos, so, b, halo_exchange);



    if (rank == size - 1) {
        std::cout << std::endl << "Press enter to exit..." << std::endl;
        std::cin.get();
    } else {
        std::cout << std::endl << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    return 0;
}
