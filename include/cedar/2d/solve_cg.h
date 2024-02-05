#ifndef CEDAR_2D_SOLVE_CG_H
#define CEDAR_2D_SOLVE_CG_H

#include <cedar/2d/types.h>
#include <cedar/kernels/solve_cg.h>

using real_t = cedar::real_t;
using len_t = cedar::len_t;

#include <src/2d/ftn/BMG2_SymStd_SETUP_cg_LU_gpu.f90.hpp>
#include <src/2d/ftn/BMG2_SymStd_SOLVE_cg_gpu.f90.hpp>

#include <cuda_profiler_api.h>

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_SETUP_cg_LU(real_t*, len_t*, len_t*, int*, real_t*, len_t*,len_t*,int*);
	void BMG2_SymStd_SOLVE_cg(real_t*, real_t*, len_t, len_t, real_t*, real_t*, len_t, len_t, int);
	void BMG_get_bc(int, int*);
}


namespace cedar { namespace cdr2 {

template <typename device=cedar::cpu>
class solve_cg_f90 : public kernels::solve_cg<stypes>
{
public:

	void setup(const stencil_op<five_pt> & so,
	           grid_func & ABD) override
	{
		this->setup_impl(so, ABD);
	}
	void setup(const stencil_op<nine_pt> & so,
	           grid_func & ABD) override
	{
		this->setup_impl(so, ABD);
	}

	template<class sten>
	void setup_impl(const stencil_op<sten> & so,
	                grid_func & ABD)
	{
            // std::cerr << "Running coarse grid setup on " << device::to_string() << std::endl;

		len_t nx, ny;
		int nstencil;
		len_t nabd1, nabd2;
		int ibc;

		auto & sod = const_cast<stencil_op<sten>&>(so);

		nx = so.len(0);
		ny = so.len(1);

		nstencil = stencil_ndirs<sten>::value;

		nabd1 = ABD.len(0);
		nabd2 = ABD.len(1);

		BMG_get_bc(params->per_mask(), &ibc);

                sod.template ensure<device>();
                ABD.template ensure<device>();

                if (device::is_gpu()) {
                    BMG2_SymStd_SETUP_cg_LU_gpu(sod, nx, ny, nstencil,
                                                ABD, nabd1, nabd2, ibc);
                } else {
                    BMG2_SymStd_SETUP_cg_LU(sod.data(), &nx, &ny, &nstencil,
                                            ABD.data(), &nabd1, &nabd2, &ibc);
                }

                // ABD.ensure_cpu();

                // std::cerr << "Coarse grid factorization: " << std::endl << ABD.to_buffer() << std::endl
	}

        void run(grid_func & x,
                 const grid_func & b,
                 const grid_func & ABD,
                 array<real_t, 1>& bbd) override
        {
            // std::cerr << "Running coarse grid solve on " << device::to_string() << std::endl;

            int ibc;

            grid_func & bd = const_cast<grid_func&>(b);
            grid_func & abd_data = const_cast<grid_func&>(ABD);

            BMG_get_bc(params->per_mask(), &ibc);

            x.template ensure<device>();
            bd.template ensure<device>();
            abd_data.template ensure<device>();
            bbd.template ensure<device>();

            auto nabd1 = ABD.len(0);
            auto nabd2 = ABD.len(1);

            if (device::is_gpu()) {
                BMG2_SymStd_SOLVE_cg_gpu(x, bd, x.len(0), x.len(1),
                                         abd_data, bbd, nabd1, nabd2, ibc);
            } else {
                BMG2_SymStd_SOLVE_cg(x.data(), bd.data(), x.len(0), x.len(1),
                                     abd_data.data(), bbd.data(), ABD.len(0), ABD.len(1), ibc);
            }

            // auto rhs_buf = std::static_pointer_cast<ftl::CUDABaseBuffer<real_t, len_t>>(bd.to_flat_buffer().get_dev_impl())->get_device_pointer();
            // const std::size_t numel = bd.to_flat_buffer().get_numel();
            // std::vector<real_t> rhs(numel);
            // cudaMemcpy(rhs.data(), rhs_buf, numel * sizeof(real_t), cudaMemcpyDeviceToHost);

            // std::cerr << "GPU RHS for coarse solve (in solve_cg.h)" << std::endl;
            // for (std::size_t i = 0; i < rhs.size(); ++i) {
            //     std::cerr << rhs[i] << " ";
            // }
            // std::cerr << std::endl;

            // std::cerr << "Coarse grid rhs: " << std::endl << bd.to_buffer() << std::endl;
            // std::cerr << "Coarse grid solution: " << std::endl << x.to_buffer() << std::endl;
        }
};

}}

#endif
