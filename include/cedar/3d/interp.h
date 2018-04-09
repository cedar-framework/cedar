#ifndef CEDAR_3D_INTERP_H
#define CEDAR_3D_INTERP_H

#include <cedar/3d/types.h>
#include <cedar/kernels/interp_add.h>
#include <cedar/kernels/setup_interp.h>
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

namespace cedar { namespace cdr3 {

class interp_f90 : public kernels::interp_add<cdr3::stypes>
{
	using prolong_op = cedar::cdr3::prolong_op;
	using grid_func = cedar::cdr3::grid_func;
	void run(const prolong_op & P,
	         const grid_func & coarse,
	         const grid_func & residual,
	         grid_func & fine) override;
};


class setup_interp_f90 : public kernels::setup_interp<cdr3::stypes>
{
	void run(const stencil_op<seven_pt> & fop,
	         const stencil_op<xxvii_pt> & cop,
	         prolong_op & P) override { this->run_impl(fop, cop, P); }
	void run(const stencil_op<xxvii_pt> & fop,
	         const stencil_op<xxvii_pt> & cop,
	         prolong_op & P) override { this->run_impl(fop, cop, P); }


	template<class sten>
	void run_impl(const stencil_op<sten> & fop,
	              const stencil_op<xxvii_pt> & cop,
	              prolong_op & P)
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
		array<real_t, 4> yo(fop.len(0), fop.len(1), 2, 14);
		int jpn;
		BMG_get_bc(params->per_mask(), &jpn);

		BMG3_SymStd_SETUP_interp_OI(fopd.data(), copd.data(),
		                            P.data(),
		                            fop.len(0), fop.len(1), fop.len(2),
		                            cop.len(0), cop.len(1), cop.len(2),
		                            ifd, nstencil, irelax, jpn,
		                            yo.data());
	}


	inline void store_fine_op(stencil_op<seven_pt> & fop,
	                          prolong_op & P)
	{
		P.fine_op_seven = &fop;
		P.fine_is_seven = true;
	}

	inline void store_fine_op(stencil_op<xxvii_pt> & fop,
	                          prolong_op & P)
	{
		P.fine_op_xxvii = &fop;
		P.fine_is_seven = false;
	}

};

}}

#endif
