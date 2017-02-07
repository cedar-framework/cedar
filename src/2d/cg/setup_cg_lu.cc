#include "boxmg/2d/cg/setup_cg_lu.h"

extern "C" {
	using namespace boxmg;
	void BMG2_SymStd_SETUP_cg_LU(real_t*, len_t*, len_t*, int*, real_t*, len_t*,len_t*,int*);
}

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	using namespace boxmg::bmg2d;
	void setup_cg_lu(const stencil_op & so,
	                 grid_func & ABD)
	{
		len_t nx, ny;
		int nstencil;
		len_t nabd1, nabd2;
		int ibc = 0;

		const grid_stencil & so_sten = so.stencil();
		stencil_op & sod = const_cast<stencil_op&>(so);

		nx = so_sten.len(0);
		ny = so_sten.len(1);

		if (so_sten.five_pt()) nstencil = 3;
		else nstencil = 5;

		nabd1 = ABD.len(0);
		nabd2 = ABD.len(1);

		BMG2_SymStd_SETUP_cg_LU(sod.data(), &nx, &ny, &nstencil,
		                        ABD.data(), &nabd1, &nabd2, &ibc);
	}
}

}}}
