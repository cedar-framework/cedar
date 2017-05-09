#include <cedar/3d/inter/setup_interp.h>
#include <cedar/2d/ftn/BMG_parameters_c.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SETUP_interp_OI(int kgf, int kgc, real_t *so, real_t *soc,
	                                 real_t *ci, len_t iif, len_t jjf, len_t kkf,
	                                 len_t iic, len_t jjc, len_t kkc,
	                                 int nog, int ifd, int nstncl, int irelax, int jpn,
	                                 real_t *yo);
	void BMG_get_bc(int, int*);
}


namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	using namespace cedar::cdr3;

	void setup_interp(const kernel_params & params,
	                  int kf, int kc, int nog,
	                  const stencil_op & fop,
	                  const stencil_op & cop,
	                  inter::prolong_op &P)
	{
		int nstencil;
		int irelax = BMG_RELAX_SYM;
		int ifd;

		const grid_stencil &fsten = fop.stencil();
		const grid_stencil &csten = cop.stencil();
		stencil_op &fopd = const_cast<stencil_op&>(fop);
		stencil_op &copd = const_cast<stencil_op&>(cop);

		P.fine_op = &fopd;

		if (fsten.five_pt()) {
			ifd = 1;
			nstencil = 4;
		} else {
			ifd = 0;
			nstencil = 14;
		}

		// TODO: preallocate this?
		array<len_t, real_t, 4> yo(fsten.len(0), fsten.len(1), 2, 14);
		int jpn;
		BMG_get_bc(params.per_mask(), &jpn);

		BMG3_SymStd_SETUP_interp_OI(kf, kc,
		                            fopd.data(), copd.data(),
		                            P.data(),
		                            fsten.len(0), fsten.len(1), fsten.len(2),
		                            csten.len(0), csten.len(1), csten.len(2),
		                            nog, ifd, nstencil, irelax, jpn,
		                            yo.data());
	}
}

}}}
