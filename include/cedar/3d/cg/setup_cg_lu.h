#ifndef CEDAR_3D_SETUP_CG_LU_H
#define CEDAR_3D_SETUP_CG_LU_H

#include "cedar/2d/ftn/BMG_parameters_c.h"
#include <cedar/kernel_params.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/grid_func.h>


extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SETUP_cg_LU(real_t *so, len_t ii, len_t jj, len_t kk, int NStncl,
	                             real_t *abd, len_t nabd1, len_t nabd2, int ibc);
	void BMG_get_bc(int, int*);
}


namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	template<class sten>
	void setup_cg_lu(const kernel_params & params,
	                 const stencil_op<sten> &so,
	                 grid_func &ABD)
	{
		int nstencil, ibc;
		auto & sod = const_cast<stencil_op<sten>&>(so);

		nstencil = stencil_ndirs<sten>::value;

		BMG_get_bc(params.per_mask(), &ibc);

		BMG3_SymStd_SETUP_cg_LU(sod.data(), so.len(0), so.len(1), so.len(2),
		                        nstencil, ABD.data(), ABD.len(0), ABD.len(1), ibc);
	}
}

}}}

#endif
