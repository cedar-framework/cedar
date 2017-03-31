#include <cedar/2d/ftn/BMG_parameters_c.h>

#include <cedar/3d/inter/interp.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_interp_add(real_t *q, real_t *qc,
	                            real_t *so, real_t *res, real_t *ci,
	                            len_t iic, len_t jjc, len_t kkc,
	                            len_t iif, len_t jjf, len_t kkf,
	                            int NStncl, int jpn);
}

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	void fortran_interp(const inter::prolong_op & P,
	                    const grid_func & coarse,
	                    const grid_func & residual,
	                    grid_func & fine)
	{
		using namespace cedar::cdr3;

		int nstencil, ibc;

		inter::prolong_op & Pd = const_cast<inter::prolong_op&>(P);
		grid_func & coarsed = const_cast<grid_func&>(coarse);
		grid_func & res = const_cast<grid_func&>(residual);

		if (Pd.stencil().five_pt()) {
			nstencil = 4;
		} else {
			nstencil = 14;
		}

		ibc = BMG_BCs_definite;

		BMG3_SymStd_interp_add(fine.data(), coarsed.data(), Pd.fine_op->data(), res.data(), Pd.data(),
		                       coarsed.len(0), coarsed.len(1), coarsed.len(2),
		                       fine.len(0), fine.len(1), fine.len(2),
		                       nstencil, ibc);
	}
}

}}}
