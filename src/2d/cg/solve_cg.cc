#include "boxmg/2d/ftn/BMG_parameters_c.h"

#include "boxmg/2d/cg/solve_cg.h"


extern "C" {
	using namespace boxmg;
	void BMG2_SymStd_SOLVE_cg(real_t*, real_t*, len_t, len_t, real_t*, real_t*, len_t, len_t, int);
}


namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void fortran_solve_cg(grid_func & x,
	                      const grid_func & b,
	                      const grid_func & ABD,
	                      real_t *bbd)
	{
		using namespace boxmg::bmg2d;

		int ibc;

		grid_func & bd = const_cast<grid_func&>(b);
		grid_func & abd_data = const_cast<grid_func&>(ABD);


		ibc = BMG_BCs_definite;

		BMG2_SymStd_SOLVE_cg(x.data(), bd.data(), x.len(0), x.len(1),
		                    abd_data.data(), &bbd[0], ABD.len(0), ABD.len(1), ibc);
	}
}

}}}
