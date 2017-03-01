#include "boxmg/2d/ftn/BMG_parameters_c.h"
#include "boxmg/2d/inter/galerkin_prod.h"


extern "C" {
	using namespace boxmg;
	void BMG2_SymStd_SETUP_ITLI_ex(int kf, int kc, real_t *so, real_t *soc, real_t *ci,
	                               len_t iif, len_t jjf, len_t iic, len_t jjc, int nog,
	                               int ifd, int nstncl, int ipn);
}


namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	using namespace boxmg::bmg2d;

	void galerkin_prod(int kf, int kc, int nog,
	                   const inter::prolong_op & P,
	                   const stencil_op & fop,
	                   stencil_op & cop)
	{
		using namespace boxmg::bmg2d;
		len_t iif, jjf, iic, jjc;
		int nstencil, ipn, ifd;

		const grid_stencil &fsten = fop.stencil();
		grid_stencil & csten = cop.stencil();
		stencil_op &fopd = const_cast<stencil_op&>(fop);
		inter::prolong_op &Pd = const_cast<inter::prolong_op&>(P);

		iif = fsten.len(0);
		jjf = fsten.len(1);
		iic = csten.len(0);
		jjc = csten.len(1);

		if (fsten.five_pt()) {
			ifd = 1;
			nstencil = 5;
		} else {
			ifd = 0;
			nstencil = 9;
		}

		ipn = BMG_BCs_definite;

		BMG2_SymStd_SETUP_ITLI_ex(kf, kc, fopd.data(), cop.data(), Pd.data(),
		                          iif, jjf, iic, jjc, nog, ifd, nstencil, ipn);

	}
}

}}}
