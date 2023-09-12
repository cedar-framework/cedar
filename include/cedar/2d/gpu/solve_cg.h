#ifndef CEDAR_2D_GPU_SOLVE_CG_H
#define CEDAR_2D_GPU_SOLVE_CG_H

#include <cedar/2d/gpu/types.h>
#include <cedar/kernels/solve_cg.h>
#include <cedar/2d/ftn/BMG_parameters_c.h>

using real_t = cedar::real_t;
using len_t = cedar::len_t;

#include <src/2d/ftn/mpi/BMG2_SymStd_GS_solve_cg.f90.hpp>

namespace cedar::cdr2::gpu::mpi {

    class solve_cg_f90 : public kernels::solve_cg<stypes>
    {
    private:
        int nstencil;
        int ifd;

    public:

	void setup(const stencil_op<five_pt> & so,
	           grid_func & ABD) override
	{
            nstencil = 3;
            ifd = 1;
            this->setup_impl(so, ABD);
	}
	void setup(const stencil_op<nine_pt> & so,
	           grid_func & ABD) override
	{
            nstencil = 5;
            ifd = 0;
            this->setup_impl(so, ABD);
	}

	template<class sten>
	void setup_impl(const stencil_op<sten> & so,
	                grid_func & ABD) {

            auto sob = so.to_buffer();
            // std::cerr << " == Setup coarse solve ==" << std::endl;
            // std::cerr << "Operator: " << sob << std::endl;
            // len_t nx, ny;
            // int nstencil;
            // len_t nabd1, nabd2;
            // int ibc;

            // auto & sod = const_cast<stencil_op<sten>&>(so);

            // len_t nx = so.len(0);
            // len_t ny = so.len(1);

            // int nstencil = stencil_ndirs<sten>::value;

            // nabd1 = ABD.len(0);
            // nabd2 = ABD.len(1);

            // BMG_get_bc(params->per_mask(), &ibc);

            // // BMG2_SymStd_SETUP_cg_LU(sod.data(), &nx, &ny, &nstencil,
            // //                         ABD.data(), &nabd1, &nabd2, &ibc);
            // MPI_BMG2_SymStd_SETUP_recip<ftl::device::GPU>(
            //     so, abd, nx, ny, nstencil);
            ABD.ensure_gpu();

            copy_so_to_abd<ftl::device::GPU>(so, ABD, so.len(0), so.len(1), so.len(2));

            auto abdb = ABD.to_buffer();
            // std::cerr << "ABD: " << std::endl << abdb << std::endl;
            // std::cerr << " ========================" << std::endl;
	}

	void run(grid_func & x,
	         const grid_func & b,
	         const grid_func & ABD,
	         real_t * bbd) override {
            int ibc;

            grid_func & bd = const_cast<grid_func&>(b);
            grid_func & abd_data = const_cast<grid_func&>(ABD);
            grid_topo & topo = x.grid();

            int k = topo.level()+1;
            int kf = topo.nlevel();
            int irelax_sym = 1;
            int updown = BMG_DOWN;
            void* halof = services->fortran_handle<halo_exchange>();

            auto bb = bd.to_buffer();
            auto xb = x.to_buffer();
            auto abdb = abd_data.to_buffer();

            std::cerr << " == Coarse-grid solve == " << std::endl;
            // std::cerr << "Operator: " << std::endl << abdb << std::endl;
            std::cerr << "RHS: " << std::endl << bb << std::endl;
            std::cerr << "Pre-solve solution:" << std::endl << xb << std::endl;

            for (int i = 0; i < 1; ++i) {
            MPI_BMG2_SymStd_GS_solve_cg<ftl::device::GPU>(
                k, abd_data, bd, x, ABD.len(0), ABD.len(1), kf, ifd, nstencil,
                irelax_sym, updown, topo.is(0), topo.is(1), halof);
            }

            std::cerr << "Post-solve solution:" << std::endl << xb << std::endl;
            std::cerr << " ======================= " << std::endl;

                // BMG_get_bc(params->per_mask(), &ibc);

                // BMG2_SymStd_SOLVE_cg(x.data(), bd.data(), x.len(0), x.len(1),
                //                      abd_data.data(), &bbd[0], ABD.len(0), ABD.len(1), ibc);

                // services->get<halo_exchange>().run(abd_data);
                }
        };

    }

#endif
