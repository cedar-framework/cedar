#include "cedar/2d/cg/setup_cg_lu.h"

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_SETUP_cg_LU(real_t*, len_t*, len_t*, int*, real_t*, len_t*,len_t*,int*);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	using namespace cedar::cdr2;
	void setup_cg_lu(const kernel_params & params, const stencil_op & so,
	                 grid_func & ABD)
	{
		len_t nx, ny;
		int nstencil;
		len_t nabd1, nabd2;
		int ibc;

		const grid_stencil & so_sten = so.stencil();
		stencil_op & sod = const_cast<stencil_op&>(so);

		nx = so_sten.len(0);
		ny = so_sten.len(1);

		if (so_sten.five_pt()) nstencil = 3;
		else nstencil = 5;

		nabd1 = ABD.len(0);
		nabd2 = ABD.len(1);

		BMG_get_bc(params.per_mask(), &ibc);

		BMG2_SymStd_SETUP_cg_LU(sod.data(), &nx, &ny, &nstencil,
		                        ABD.data(), &nabd1, &nabd2, &ibc);
	}
}

}}}
