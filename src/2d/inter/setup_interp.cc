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

	template<>
	void setup_interp(const kernel_params & params,
	                  int kf, int kc, int nog,
	                  const stencil_op<five_pt> &fop,
	                  const stencil_op<nine_pt> &cop,
	                  inter::prolong_op &P)
	{
		using namespace cedar::cdr2;
		int ifd, jpn;
		len_t iif, jjf, iic, jjc;
		int nstencil;
		int irelax = 0;

		auto & fopd = const_cast<stencil_op<five_pt>&>(fop);
		auto & copd = const_cast<stencil_op<nine_pt>&>(cop);

		P.fine_op_five = &fopd;
		P.fine_is_five = true;

		iif = fop.len(0);
		jjf = fop.len(1);
		iic = cop.len(0);
		jjc = cop.len(1);

		ifd = 1;
		nstencil = 3;

		BMG_get_bc(params.per_mask(), &jpn);

		BMG2_SymStd_SETUP_interp_OI(kf, kc, fopd.data(), copd.data(), P.data(), iif, jjf, iic, jjc,
		                            nog, ifd, nstencil, jpn, irelax);
	}


	template<>
	void setup_interp(const kernel_params & params,
	                  int kf, int kc, int nog,
	                  const stencil_op<nine_pt> &fop,
	                  const stencil_op<nine_pt> &cop,
	                  inter::prolong_op &P)
	{
		using namespace cedar::cdr2;
		int ifd, jpn;
		len_t iif, jjf, iic, jjc;
		int nstencil;
		int irelax = 0;

		auto & fopd = const_cast<stencil_op<nine_pt>&>(fop);
		auto & copd = const_cast<stencil_op<nine_pt>&>(cop);

		P.fine_op_nine = &fopd;
		P.fine_is_five = false;

		iif = fop.len(0);
		jjf = fop.len(1);
		iic = cop.len(0);
		jjc = cop.len(1);

		ifd = 0;
		nstencil = 5;

		BMG_get_bc(params.per_mask(), &jpn);

		BMG2_SymStd_SETUP_interp_OI(kf, kc, fopd.data(), copd.data(), P.data(), iif, jjf, iic, jjc,
		                            nog, ifd, nstencil, jpn, irelax);
	}
}

}}}
