#ifndef CEDAR_2D_COARSEN_H
#define CEDAR_2D_COARSEN_H

#include <cedar/kernels/coarsen_op.h>
#include <cedar/2d/types.h>
#include <cedar/2d/ftn/BMG_parameters_c.h>

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_SETUP_ITLI_ex(real_t *so, real_t *soc, real_t *ci,
	                               len_t iif, len_t jjf, len_t iic, len_t jjc,
	                               int ifd, int nstncl, int ipn);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr2 {

class galerkin : public kernels::coarsen_op<stypes>
{
	void run(const prolong_op & P,
	         const stencil_op<five_pt> & fop,
	         stencil_op<nine_pt> & cop) override
	{
		this->run_impl(P, fop, cop);
	}
	void run(const prolong_op & P,
	         const stencil_op<nine_pt> & fop,
	         stencil_op<nine_pt> & cop) override
	{
		this->run_impl(P, fop, cop);
	}


	template<class sten>
	void run_impl(const prolong_op & P,
	              const stencil_op<sten> & fop,
	              stencil_op<nine_pt> & cop)
	{
		len_t iif, jjf, iic, jjc;
		int nstencil, ipn, ifd;

		auto & fopd = const_cast<stencil_op<sten>&>(fop);
		prolong_op &Pd = const_cast<prolong_op&>(P);

		iif = fop.len(0);
		jjf = fop.len(1);
		iic = cop.len(0);
		jjc = cop.len(1);

		// this may have to be 5 and 9 instead of 3 and 5
		nstencil = stencil_ndirs<sten>::value;
		if (std::is_same<sten, five_pt>::value)
			ifd = 1;
		else
			ifd = 0;

		BMG_get_bc(params->per_mask(), &ipn);

		BMG2_SymStd_SETUP_ITLI_ex(fopd.data(), cop.data(), Pd.data(),
		                          iif, jjf, iic, jjc, ifd, nstencil, ipn);
	}


};


}}

#endif
