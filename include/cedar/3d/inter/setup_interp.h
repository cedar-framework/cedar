#ifndef CEDAR_3D_INTER_SETUP_INTERP_H
#define CEDAR_3D_INTER_SETUP_INTERP_H

#include <type_traits>

#include <cedar/kernel_params.h>
#include <cedar/3d/grid_func.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/inter/prolong_op.h>

#include <cedar/2d/ftn/BMG_parameters_c.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SETUP_interp_OI(real_t *so, real_t *soc,
	                                 real_t *ci, len_t iif, len_t jjf, len_t kkf,
	                                 len_t iic, len_t jjc, len_t kkc,
	                                 int ifd, int nstncl, int irelax, int jpn,
	                                 real_t *yo);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	template<class sten>
		void store_fine_op(stencil_op<sten> & fop,
		                   inter::prolong_op & P);
	template<>
		void store_fine_op(stencil_op<seven_pt> & fop,
		                   inter::prolong_op & P)
	{
		P.fine_op_seven = &fop;
		P.fine_is_seven = true;
	}

	template<>
		void store_fine_op(stencil_op<xxvii_pt> & fop,
		                   inter::prolong_op & P)
	{
		P.fine_op_xxvii = &fop;
		P.fine_is_seven = false;
	}

	template<class sten>
	void setup_interp(const kernel_params & params,
	                  const stencil_op<sten> & fop,
	                  const stencil_op<xxvii_pt> & cop,
	                  inter::prolong_op & P)
	{
		int nstencil;
		int irelax = BMG_RELAX_SYM;
		int ifd;

		auto &fopd = const_cast<stencil_op<sten>&>(fop);
		auto &copd = const_cast<stencil_op<xxvii_pt>&>(cop);

		store_fine_op(fopd, P);

		nstencil = stencil_ndirs<sten>::value;

		if (std::is_same<sten, seven_pt>::value)
			ifd = 1;
		else
			ifd = 0;

		// TODO: preallocate this?
		array<len_t, real_t, 4> yo(fop.len(0), fop.len(1), 2, 14);
		int jpn;
		BMG_get_bc(params.per_mask(), &jpn);

		BMG3_SymStd_SETUP_interp_OI(fopd.data(), copd.data(),
		                            P.data(),
		                            fop.len(0), fop.len(1), fop.len(2),
		                            cop.len(0), cop.len(1), cop.len(2),
		                            ifd, nstencil, irelax, jpn,
		                            yo.data());
	}
}

}}}

#endif
