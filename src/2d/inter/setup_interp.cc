#include "cedar/2d/inter/setup_interp.h"

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_SETUP_interp_OI(int kf, int kc, real_t *so, real_t *soc, real_t *ci,
	                                 len_t iif, len_t jjf, len_t iic, len_t jjc,
	                                 int nog, int ifd, int nstncl, int jpn, int irelax);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	using namespace cedar::cdr2;

	void setup_interp(const kernel_params & params,
	                  int kf, int kc, int nog,
	                  const stencil_op &fop,
	                  const stencil_op &cop,
	                  inter::prolong_op &P)
	{
		using namespace cedar::cdr2;
		int ifd, jpn;
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

		BMG_get_bc(params.per_mask(), &jpn);

		// std::cout << "fine stencil: " << fopd.data() << std::endl;
		// std::cout << "coarse stencil: " << copd.data() << std::endl;
		// std::cout << "interpolation op: " << P.data() << std::endl;
		BMG2_SymStd_SETUP_interp_OI(kf, kc, fopd.data(), copd.data(), P.data(), iif, jjf, iic, jjc,
		                            nog, ifd, nstencil, jpn, irelax);
	}
}

}}}
