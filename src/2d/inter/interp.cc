#include "boxmg/2d/ftn/BMG_parameters_c.h"
#include "boxmg/2d/inter/interp.h"

extern "C" {
	using namespace boxmg;
	void BMG2_SymStd_interp_add(real_t*, real_t*, real_t*,real_t*,real_t*,
	                            len_t, len_t, len_t, len_t, int, int);
}

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void fortran_interp(const inter::prolong_op & P,
	                    const grid_func & coarse,
	                    const grid_func & residual,
	                    grid_func & fine)
	{
		using namespace boxmg::bmg2d;

		int nstencil, ibc;

		inter::prolong_op & Pd = const_cast<inter::prolong_op&>(P);
		grid_func & coarsed = const_cast<grid_func&>(coarse);
		grid_func & res = const_cast<grid_func&>(residual);

		if (Pd.stencil().five_pt()) {
			nstencil = 3;
		} else {
			nstencil = 5;
		}

		ibc = BMG_BCs_definite;

		BMG2_SymStd_interp_add(fine.data(), coarsed.data(), res.data(), Pd.fine_op->data(), Pd.data(),
		                       coarsed.len(0), coarsed.len(1), fine.len(0), fine.len(1),
		                       nstencil, ibc);
	}
}


}}}
