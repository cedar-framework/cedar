#ifndef VCYCLE_HPP_INC_
#define VCYCLE_HPP_INC_

#include <ftl/KernelRegistry.hpp>
#include <ftl/Buffer.hpp>
#include <ftl/Runtime.hpp>
#include <ftl/String.hpp>
#include <ftl/Device.hpp>
#include <ftl/Cedar.hpp>

#include <random>

#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_recip.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_interp_OI.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_ITLI_ex.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_residual.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_restrict.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_interp_add.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_relax_GS.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_fine_stencil.f90.hpp>

#include <cedar/2d/mpi/kernel_manager.h>
#include <cedar/2d/ftn/mpi/BMG_workspace_c.h>

/*
  #include <cedar/services/halo_exchange.h>
  #include <cedar/mpi/grid_topo.h>
*/

using grid_dim = std::tuple<int, int>;

template <typename T, typename Device=ftl::device::GPU, bool timing=true>
class VCycle {
public:
    using five_pt = cedar::cdr2::five_pt;
    std::vector<ftl::Buffer<T>> so_hierarchy;
    std::vector<ftl::Buffer<T>> sor_hierarchy;
    std::vector<ftl::Buffer<T>> interp_hierarchy;

    std::vector<cedar::topo_ptr> topos;

private:

    cedar::cdr2::mpi::stencil_op<five_pt>* so_fine_cedar;

    /* Pre-allocate space for rhs and soln for v-cycle */
    std::vector<ftl::Buffer<T>> x_store;
    std::vector<ftl::Buffer<T>> b_store;
    std::vector<ftl::Buffer<T>> r_store;

    MPI_Comm communicator;
    cedar::cdr2::mpi::kman_ptr kman;

    ftl::Buffer<T> u_fine;

    int jpn;
    int irelax;

    int smoothing;

    std::size_t num_levels;

    void create_smoother_coeff(ftl::Buffer<T>& sor, ftl::Buffer<T>& so) {
        const auto shape = so.get_shape();
        MPI_BMG2_SymStd_SETUP_recip<Device>(so, sor, shape[0], shape[1], shape[2]);
        std::cout << "MPI_BMG2_SymStd_SETUP_recip(so, sor, " << shape[0] << ", " << shape[1] << ", " << shape[2] << ");" << std::endl;
        //std::cout << sor << std::endl;
    }

    void push_back_smoother_coeff(ftl::Buffer<T>& so) {
        const auto shape = so.get_shape();
        ftl::Buffer<T> sor({shape[0], shape[1], 3}, ftl::Buffer_Both);
        create_smoother_coeff(sor, so);
        sor_hierarchy.push_back(sor);
    }

public:

    std::pair<int, int> get_coarse_size(int iif, int jjf) {
        const int iic = (iif - 1)/2 + 2;
        const int jjc = (jjf - 1)/2 + 2;
        return std::make_pair(iic, jjc);
    }

    void add_level(std::size_t level) {
        assert(level >= 0 && level < num_levels -1);

        /* Get dimensions of coarsest level */
        const int global_i = topos[level]->nglobal(0);
        const int global_j = topos[level]->nglobal(1);
        const int local_i = topos[level]->nlocal(0);
        const int local_j = topos[level]->nlocal(1);

        const int global_ic = topos[level + 1]->nglobal(0);
        const int global_jc = topos[level + 1]->nglobal(1);
        const int local_ic = topos[level + 1]->nlocal(0);
        const int local_jc = topos[level + 1]->nlocal(1);

        /* Fine grid topology */
        cedar::topo_ptr topo = topos[level];

        ftl::Buffer<T>& so = so_hierarchy.back();
        const int nstncl = so.get_shape()[2];
        const int ifd = int(nstncl == 3);

        /* First, create the interpolation operator */
        ftl::Buffer<T> soc({local_ic + 1, local_jc + 1, 5}, ftl::Buffer_Both);
        ftl::Buffer<T> ci({local_ic, local_jc, 8}, ftl::Buffer_Both);

        int kc = topos[level]->level();
        int kf = kc + 1;
        int nog = topo->nlevel();

        cedar::grid_topo::igrd_t igrd = topo->get_igrd();
        ftl::Buffer<uint32_t> igrd_buf(igrd->data(), {nog, NBMG_pIGRD});
        void* halof = halo_exchange();

        MPI_BMG2_SymStd_SETUP_interp_OI<Device>(
            kf, kc, so, ci, local_i, local_j, local_ic, local_jc,
            nog, nog, igrd_buf, ifd, nstncl, jpn, halof);

        std::cout << "MPI_BMG2_SymStd_SETUP_interp_OI(" << kf << ", " << kc << ", fopd, P, " << local_i << ", " << local_j << ", " <<
            local_ic << ", " << local_jc << ", " << nog << ", " << nog << ", igrd, " << ifd << ", " << nstncl << ", " << jpn <<
            ", halof)" << std::endl;

        //std::cout << ci << std::endl;

        /* With the interpolation operator, create the coarse-grid operator */
        MPI_BMG2_SymStd_SETUP_ITLI_ex<Device>(
            kf, kc, so, soc, ci, local_i, local_j, local_ic, local_jc,
            topo->is(0), topo->is(1), nog, ifd, nstncl, halof);

        /* Create the smoother coefficients */
        push_back_smoother_coeff(soc);

        so_hierarchy.push_back(soc);
        interp_hierarchy.push_back(ci);

        const std::vector<int> vec_shape = {local_ic, local_jc};
        x_store.emplace_back(vec_shape, ftl::Buffer_Both);
        b_store.emplace_back(vec_shape, ftl::Buffer_Both);
        r_store.emplace_back(vec_shape, ftl::Buffer_Both);
    }

    void setup_halo() {
        auto& sman = kman->services();
	auto& halo_service = sman.get<cedar::cdr2::mpi::halo_exchange>();
	halo_service.setup(topos);
    }

    void setup_space() {
        for (std::size_t i = 0; i < num_levels - 1; ++i) {
            auto& fine_grid = topos[i];

            if (i == 0) {
                fine_grid->comm = communicator;
                fine_grid->grow(num_levels);
            }

            const int kc = num_levels - i - 2;

            auto coarse_grid = std::make_shared<cedar::grid_topo>(
                fine_grid->get_igrd(), kc, num_levels);
            coarse_grid->comm = fine_grid->comm;

            len_t NLxg = fine_grid->nlocal(0) - 2;
            len_t NLyg = fine_grid->nlocal(1) - 2;
            len_t NGxg = (fine_grid->nglobal(0) - 1) / 2 + 2;
            len_t NGyg = (fine_grid->nglobal(1) - 1) / 2 + 2;

            coarse_grid->nglobal(0) = NGxg;
            coarse_grid->nglobal(1) = NGyg;

            if ((fine_grid->is(0) % 2) == 1) {
                coarse_grid->is(0) = (fine_grid->is(0) + 1) / 2;
                NLxg = (NLxg + 1) / 2;
            } else {
                coarse_grid->is(0) = fine_grid->is(0)/2 + 1;
                if (NLxg % 2 == 1) {
                    NLxg = (NLxg-1)/2;
                } else {
                    NLxg = (NLxg+1)/2;
                }
            }

            if (fine_grid->is(1) % 2 == 1) {
                coarse_grid->is(1) = (fine_grid->is(1)+1) / 2;
                NLyg = (NLyg+1) / 2;
            } else {
                coarse_grid->is(1) = fine_grid->is(1) / 2 + 1;
                if (NLyg % 2 == 1) {
                    NLyg = (NLyg - 1) / 2;
                } else {
                    NLyg = (NLyg+1)/2;
                }
            }

            coarse_grid->nlocal(0) = NLxg + 2;
            coarse_grid->nlocal(1) = NLyg + 2;

            coarse_grid->nproc(0) = fine_grid->nproc(0);
            coarse_grid->nproc(1) = fine_grid->nproc(1);
            coarse_grid->coord(0) = fine_grid->coord(0);
            coarse_grid->coord(1) = fine_grid->coord(1);

            topos.push_back(coarse_grid);
	}

	setup_halo();

        //auto service = kman->services().fortran_handle<cedar::cdr2::mpi::halo_exchange>();
        //service->exchange_sten(1, so_hierarchy[0].data());
    }

    void initialize() {
        int i, j;
        i = topos[0]->nglobal(0);
        j = topos[0]->nglobal(1);

        /* Count how many levels we will have overall for the full hierarchy */
        while (i > 6 && j > 6) {
            std::tie(i, j) = get_coarse_size(i, j);
            num_levels ++;
        }

        /* Setup cedar grid topologies */
        setup_space();

        for (int i = 0; i < num_levels - 1; ++i) {
            add_level(i);
        }

        for (int i = 0; i < num_levels; ++i) {
            std::cout << " == " << i << " == " << std::endl;
            std::cout << " SO " << std::endl;
            std::cout << so_hierarchy[i] << std::endl;
            std::cout << " SOR " << std::endl;
            std::cout << sor_hierarchy[i] << std::endl;

            if (i < num_levels - 1) {
                std::cout << " CI " << std::endl;
                std::cout << interp_hierarchy[i] << std::endl;
            }
        }
    }

    void* halo_exchange() {
        auto ptr = kman->services().fortran_handle<cedar::cdr2::mpi::halo_exchange>();
        return static_cast<void*>(ptr);
    }

    template <typename CycleDevice=Device>
    void cycle(ftl::Buffer<T>& b, ftl::Buffer<T>& x, int level=0) {
        const bool coarsest = (level == so_hierarchy.size() - 1);

        ftl::Buffer<T>& so = so_hierarchy[level];

        auto topo = topos[level];
        int global_i = topo->nglobal(0);
        int global_j = topo->nglobal(1);
        int local_i = topo->nlocal(0);
        int local_j = topo->nlocal(1);

        int nstncl = so.get_shape()[2];
        int ifd = int(nstncl == 3);

        int kf = topo->level() + 1;
        int kc = kf - 1;
        int nog = kf;
        int nlevel = topo->nlevel();

        // std::cout << "Level " << level << std::endl;
        // std::cout << "global_i: " << global_i << std::endl;
        // std::cout << "global_j: " << global_j << std::endl;
        // std::cout << "local_i: " << local_i << std::endl;
        // std::cout << "local_j: " << local_j << std::endl;
        // std::cout << "nstncl: " << nstncl << std::endl;
        // std::cout << "ifd: " << ifd << std::endl;
        // std::cout << "kf: " << kf << std::endl;
        // std::cout << "kc: " << kc << std::endl;
        // std::cout << "nog: " << nog << std::endl;
        // std::cout << "nlevel: " << nlevel << std::endl;
        // std::cout << "=========" << std::endl << std::endl;

        int irelax_sym = 0;
        int mpicomm = 0;

        void* halof = halo_exchange();

        const std::string level_str = "(" + std::to_string(level) + ") ";

        ftl::Buffer<T>& sor = sor_hierarchy[level];

        /* Pre-smooth */
        for (std::size_t i = 0; i < smoothing; ++i) {
            MPI_BMG2_SymStd_relax_GS
                <CycleDevice>(kf, so, b, x, sor, local_i, local_j, nlevel,
                              ifd, nstncl, 3, BMG_RELAX_SYM, BMG_DOWN, jpn,
                              halof);
        }

        std::cout << "Pre-Smoothing " << std::endl << x << std::endl;

        /* Recurse */
        if (!coarsest) {
            /* Coarse-level variables */
            ftl::Buffer<T>& ci = interp_hierarchy[level];

            auto topo_c = topos[level + 1];
            int local_ic = topo_c->nlocal(0);
            int local_jc = topo_c->nlocal(1);

            /* Compute residual */
            ftl::Buffer<T>& r = r_store[level];
            ftl::Buffer<T>& r_H = b_store[level + 1];
            ftl::Buffer<T>& x_H = x_store[level + 1];

            MPI_BMG2_SymStd_residual<CycleDevice>(kf, nlevel, nog, so, b, x, r, local_i, local_j, ifd, nstncl, 0, 0, 0);

            std::cout << "Residual " << std::endl << r << std::endl;

            /* Restrict residual */
            MPI_BMG2_SymStd_restrict<CycleDevice>(kf, 1, nog, r, r_H, ci, local_i, local_j, local_ic, local_jc, topo->is(0), topo->is(1));

            std::cout << "Coarse residual " << std::endl << r_H << std::endl;

            /* Recurse */
            // if (local_ic <= 64 || local_jc <= 128) {
            //     x_H.template zero_on<ftl::device::CPU>();
            //     cycle<ftl::device::CPU>(r_H, x_H, level + 1);
            // } else {
            //x_H.template zero_on<CycleDevice>();
            x_H.zero();
            cycle<CycleDevice>(r_H, x_H, level + 1);
            // }

            /* Interpolate error */
            MPI_BMG2_SymStd_interp_add<CycleDevice>(
                kc, kf, nog, x, x_H, r, so, nstncl, ci,
                local_ic, local_jc, local_i, local_j,
                topo->is(0), topo->is(1), halof);

            std::cout << "Interpolated error " << std::endl << x << std::endl;
        } else {
            /* TODO: Perform an actual coarse-level solve, for now we do an absurd
               number of smoothing steps.*/

            for (std::size_t i = 0; i < 10; ++i) {
                MPI_BMG2_SymStd_relax_GS
                    <CycleDevice>(kf, so, b, x, sor, local_i, local_j, 1,
                                  ifd, nstncl, 3, BMG_RELAX_SYM, BMG_UP, jpn,
                                  halof);
            }
        }

        /* Post-smooth */
        for (std::size_t i = 0; i < smoothing; ++i) {
            MPI_BMG2_SymStd_relax_GS
                <CycleDevice>(kf, so, b, x, sor, local_i, local_j, 1,
                              ifd, nstncl, 3, BMG_RELAX_SYM, BMG_UP, jpn,
                              halof);
        }

        std::cout << "Post-Smoothing " << std::endl << x << std::endl;
    }

    ftl::Buffer<T>& gen_x0(int seed=0) {
        std::default_random_engine rng(seed);
        std::normal_distribution<real_t> normal;
        ftl::Buffer<T>& x = x_store[0];

        const auto shape = x.get_shape();
        const int iif = shape[0];
        const int jjf = shape[1];

        /* Generate a random rhs */
        for (int i = 0; i < iif; ++i) {
            for (int j = 0; j < jjf; ++j) {

                if (i != 0 && i != iif-1 &&
                    j != 0 && j != jjf-1) {
                    x(i, j) = normal(rng);
                }
            }
        }

        // if (Device::is_device()) {
        //     x.host_to_dev();
        // }

        return x;
    }

    ftl::Buffer<T>& solve(int seed=0) {
        auto x = gen_x0(seed);
        cycle(b_store[0], x, 0);
        return x;
    }

    void solve(ftl::Buffer<T>& x0) {
        cycle(b_store[0], x0, 0);
    }

    VCycle(cedar::config& conf, cedar::cdr2::mpi::stencil_op<five_pt>* fine_grid, ftl::Buffer<T>& u_fine, ftl::Buffer<T>& f_fine, cedar::topo_ptr grid):
        jpn(BMG_BCs_definite), irelax(0), u_fine(u_fine), smoothing(2), num_levels(1) {
        kman = cedar::cdr2::mpi::build_kernel_manager(conf);

        so_hierarchy.push_back(fine_grid);

        const auto shape = u_fine.get_shape();
        x_store.emplace_back(shape, ftl::Buffer_Both);
        b_store.push_back(f_fine);
        r_store.emplace_back(shape, ftl::Buffer_Both);

        push_back_smoother_coeff(fine_grid);

        topos.push_back(grid);

        communicator = MPI_COMM_WORLD;
    }

    template <typename Q, typename D, bool B>
    friend std::ostream& operator<<(std::ostream& out, VCycle<Q, D, B>& cycle);
};


template <typename Q, typename D, bool B>
std::ostream& operator<<(std::ostream& out, VCycle<Q, D, B>& cycle) {
    out << "boxmg solver (" << D::to_string() << ")" << std::endl;
    out << "Number of levels:\t" << cycle.so_hierarchy.size() << std::endl;
    out << "  level\t width\theight\tstencil" << std::endl;

    for (std::size_t i = 0; i < cycle.so_hierarchy.size(); ++i) {
        auto shape = cycle.so_hierarchy[i].get_shape();
        out << std::setw(5) << std::right << i << "\t"
            << std::setw(6) << std::right << shape[0] << "\t"
            << std::setw(6) << std::right << shape[1] << "\t"
            << std::setw(7) << std::right << shape[2] << std::endl;
    }
    return out;
}


#endif
