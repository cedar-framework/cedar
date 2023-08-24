#include <mpi.h>
#include <iostream>
#include <memory>
#include <math.h>

#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/mpi/gallery.h>

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
    std::vector<int32_t> shape(ND);

    for (std::size_t i = 0; i < ND; ++i) {
        const int32_t dim = array.len(i);
        shape[i] = dim;
    }

    ftl::Buffer<dtype> buf(array.data(), shape);

    return buf;
}


template <typename stencil_type>
ftl::Buffer<real_t> so_to_buffer(cedar::cdr2::mpi::stencil_op<stencil_type>& array) {
    std::vector<int32_t> shape = {array.len(0) + 1, array.len(1) + 1, array.len(2)};
    ftl::Buffer<real_t> buf(array.data(), shape);

    return buf;
}


template <typename T>
void save_file(ftl::Buffer<T>& buf, const std::string& fname) {
    std::ofstream file(fname);
    file << buf.tostr_unformatted() << std::endl;
    file.close();
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
    auto nx = ndofs[0];
    auto ny = ndofs[1];

    auto grid = util::create_topo(conf);

    /* Stencil operator */
    auto so_cedar = mpi::gallery::poisson(grid);
    ftl::Buffer<real_t> so = so_to_buffer<five_pt>(so_cedar);

    std::cout << so_cedar.len(0) << " " << so_cedar.len(1) << " " << so_cedar.len(2) << std::endl;
    std::cout << grid->nlocal(0) << " " << grid->nlocal(1) << std::endl;
    std::cout << grid->nglobal(0) << " " << grid->nglobal(1) << std::endl;

    std::cout << "so" << std::endl;
    std::cout << so.tostr_unformatted() << std::endl;
    auto* scp = so_cedar.data();
    for (std::size_t i = 0; i < 7 * 7 * 3; ++i) {
        std::cout << scp[i] << " ";
    }
    std::cout << std::endl;

    /* RHS */
    mpi::grid_func b_cedar(grid);
    set_problem(b_cedar);
    ftl::Buffer<real_t> b = to_buffer(b_cedar);

    /* Exact solution */
    mpi::grid_func exact_sol_cedar(grid);
    set_solution(exact_sol_cedar);
    ftl::Buffer<real_t> exact_sol = to_buffer(exact_sol_cedar);

    std::cout << "Set up problem" << std::endl;

    VCycle solver(conf, so, exact_sol, b, grid);

    std::cout << "Created solver object" << std::endl;

    solver.initialize();

    std::cout << "Initialized solver" << std::endl;

    std::cout << solver << std::endl;

    ftl::Buffer<real_t>& soln = solver.gen_x0(0);

    std::cout << "Generated x0" << std::endl;

    // for (int i = 0; i < 10; ++i) {
    //     solver.solve(soln);
    // }
    solver.solve(soln);

    std::cout << "Ran solver" << std::endl;

    MPI_Barrier(MPI_COMM_WORLD); // synchronize before timing solve
    MPI_Finalize();

    std::cout << "Finished" << std::endl;

    save_file(soln, "soln_" + std::to_string(rank) + ".txt");
    save_file(b, "rhs_" + std::to_string(rank) + ".txt");

    return 0;
}
