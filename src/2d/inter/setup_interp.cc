#include "boxmg/2d/inter/setup_interp.h"

extern "C" {
	using namespace boxmg;
	void BMG2_SymStd_SETUP_interp_OI(int kf, int kc, real_t *so, real_t *soc, real_t *ci,
	                                 len_t iif, len_t jjf, len_t iic, len_t jjc,
	                                 int nog, int ifd, int nstncl, int irelax);
}

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	using namespace boxmg::bmg2d;

	void setup_interp(int kf, int kc, int nog,
	                  const stencil_op &fop,
	                  const stencil_op &cop,
	                  inter::prolong_op &P)
	{
		using namespace boxmg::bmg2d;
		int ifd;
		len_t iif, jjf, iic, jjc;
		int nstencil;
		int irelax = 0;

		const grid_stencil &fsten = fop.stencil();
		const grid_stencil &csten = cop.stencil();
		stencil_op &fopd = const_cast<stencil_op&>(fop);
		stencil_op &copd = const_cast<stencil_op&>(cop);

		P.fine_op = &fopd;

		iif = fsten.len(0);
		jjf = fsten.len(1);
		iic = csten.len(0);
		jjc = csten.len(1);

		if (fsten.five_pt()) {
			ifd = 1;
			nstencil = 3;
		} else {
			ifd = 0;
			nstencil = 5;
		}

		// std::cout << "fine stencil: " << fopd.data() << std::endl;
		// std::cout << "coarse stencil: " << copd.data() << std::endl;
		// std::cout << "interpolation op: " << P.data() << std::endl;
		BMG2_SymStd_SETUP_interp_OI(kf, kc, fopd.data(), copd.data(), P.data(), iif, jjf, iic, jjc,
		                            nog, ifd, nstencil, irelax);
	}
}

}}}
