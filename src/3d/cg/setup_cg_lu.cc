#include "cedar/2d/ftn/BMG_parameters_c.h"

#include <cedar/3d/cg/setup_cg_lu.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SETUP_cg_LU(real_t *so, len_t ii, len_t jj, len_t kk, int NStncl,
	                             real_t *abd, len_t nabd1, len_t nabd2, int ibc);
	void BMG_get_bc(int, int*);
}


namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	void setup_cg_lu(const kernel_params & params, const stencil_op &so,
	                 grid_func &ABD)
	{
		int nstencil, ibc;
		auto & sod = const_cast<stencil_op&>(so);
		auto & so_sten = sod.stencil();

		ibc = BMG_BCs_definite;

		if (so_sten.five_pt()) {
			nstencil = 4;
		} else {
			nstencil = 14;
		}

		BMG_get_bc(params.per_mask(), &ibc);

		BMG3_SymStd_SETUP_cg_LU(so_sten.data(), so_sten.len(0), so_sten.len(1), so_sten.len(2),
		                        nstencil, ABD.data(), ABD.len(0), ABD.len(1), ibc);
	}
}

}}}
