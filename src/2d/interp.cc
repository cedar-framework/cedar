#include <cedar/2d/ftn/BMG_parameters_c.h>

#include <cedar/2d/interp.h>
#include <cedar/2d/types.h>

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_interp_add(real_t*, real_t*, real_t*,real_t*,real_t*,
	                            len_t, len_t, len_t, len_t, int, int);
	void BMG2_SymStd_interp_add_offload(real_t*, real_t*, real_t*,real_t*,real_t*,
	                                    len_t, len_t, len_t, len_t, int, int);
	void BMG2_SymStd_SETUP_interp_OI(real_t *so, real_t *soc, real_t *ci,
	                                 len_t iif, len_t jjf, len_t iic, len_t jjc,
	                                 int ifd, int nstncl, int jpn, int irelax);
	void BMG_get_bc(int, int*);
}


using namespace cedar;
using namespace cedar::cdr2;

interp_f90::interp_f90(bool offload)
{
	fcall = BMG2_SymStd_interp_add;
	#ifdef OFFLOAD
	if (offload)
		fcall = BMG2_SymStd_interp_add_offload;
	#endif
}

void interp_f90::run(const prolong_op & P,
                     const grid_func & coarse,
                     const grid_func & residual,
                     grid_func & fine)
{
		int nstencil, ibc;

		prolong_op & Pd = const_cast<prolong_op&>(P);
		grid_func & coarsed = const_cast<grid_func&>(coarse);
		grid_func & res = const_cast<grid_func&>(residual);

		real_t * fop_data;

		if (Pd.fine_is_five) {
			nstencil = 3;
			fop_data = Pd.fine_op_five->data();
		} else {
			nstencil = 5;
			fop_data = Pd.fine_op_nine->data();
		}

		BMG_get_bc(params->per_mask(), &ibc);

		fcall(fine.data(), coarsed.data(), res.data(), fop_data, Pd.data(),
		      coarsed.len(0), coarsed.len(1), fine.len(0), fine.len(1),
		      nstencil, ibc);
}


void setup_interp_f90::run(const stencil_op<five_pt> & fop,
                           const stencil_op<nine_pt> & cop,
                           prolong_op & P)
{
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

		BMG_get_bc(params->per_mask(), &jpn);

		BMG2_SymStd_SETUP_interp_OI(fopd.data(), copd.data(), P.data(), iif, jjf, iic, jjc,
		                            ifd, nstencil, jpn, irelax);
}

void setup_interp_f90::run(const stencil_op<nine_pt> & fop,
                           const stencil_op<nine_pt> & cop,
                           prolong_op & P)
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

		BMG_get_bc(params->per_mask(), &jpn);

		BMG2_SymStd_SETUP_interp_OI(fopd.data(), copd.data(), P.data(), iif, jjf, iic, jjc,
		                            ifd, nstencil, jpn, irelax);
}
