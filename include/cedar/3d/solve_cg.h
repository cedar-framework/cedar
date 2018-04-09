#ifndef CEDAR_3D_SOLVE_CG_H
#define CEDAR_3D_SOLVE_CG_H

#include <cedar/2d/ftn/BMG_parameters_c.h>
#include <cedar/3d/types.h>
#include <cedar/kernels/solve_cg.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SETUP_cg_LU(real_t *so, len_t ii, len_t jj, len_t kk, int NStncl,
	                             real_t *abd, len_t nabd1, len_t nabd2, int ibc);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr3 {

class solve_cg_f90 : public kernels::solve_cg<stypes>
{
public:

	void setup(const stencil_op<seven_pt> & so,
	           grid_func & ABD) override
	{
		this->setup_impl(so, ABD);
	}
	void setup(const stencil_op<xxvii_pt> & so,
	           grid_func & ABD) override
	{
		this->setup_impl(so, ABD);
	}

	template<class sten>
	void setup_impl(const stencil_op<sten> & so,
	                grid_func & ABD)
	{
		int nstencil, ibc;
		auto & sod = const_cast<stencil_op<sten>&>(so);

		nstencil = stencil_ndirs<sten>::value;

		BMG_get_bc(params->per_mask(), &ibc);

		BMG3_SymStd_SETUP_cg_LU(sod.data(), so.len(0), so.len(1), so.len(2),
		                        nstencil, ABD.data(), ABD.len(0), ABD.len(1), ibc);
	}

	void run(grid_func & x,
	         const grid_func & b,
	         const grid_func & ABD,
	         real_t * bbd) override;
};

}}

#endif
