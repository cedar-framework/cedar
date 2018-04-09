#ifndef CEDAR_3D_COARSEN_H
#define CEDAR_3D_COARSEN_H

#include <type_traits>

#include <cedar/2d/ftn/BMG_parameters_c.h>
#include <cedar/kernels/coarsen_op.h>
#include <cedar/3d/types.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SETUP_ITLI07_ex(real_t *so, real_t *soc, real_t *ci, len_t iif, len_t jjf,
	                                 len_t kkf, len_t iic, len_t jjc, len_t kkc,
	                                 int ipn);
	void BMG3_SymStd_SETUP_ITLI27_ex(real_t *so, real_t *soc, real_t *ci, len_t iif, len_t jjf,
	                                 len_t kkf, len_t iic, len_t jjc, len_t kkc, int ipn);
	void BMG_get_bc(int, int*);
}


namespace cedar { namespace cdr3 {

class galerkin : public kernels::coarsen_op<stypes>
{
	void run(const prolong_op & P,
	         const stencil_op<seven_pt> & fop,
	         stencil_op<xxvii_pt> & cop) override
	{
		this->run_impl(P, fop, cop);
	}
	void run(const prolong_op & P,
	         const stencil_op<xxvii_pt> & fop,
	         stencil_op<xxvii_pt> & cop) override
	{
		this->run_impl(P, fop, cop);
	}


	template<class sten>
	void run_impl(const prolong_op & P,
	              const stencil_op<sten> & fop,
	              stencil_op<xxvii_pt> & cop)
	{
		int ipn;
		auto &fopd = const_cast<stencil_op<sten>&>(fop);
		prolong_op &Pd = const_cast<prolong_op&>(P);

		BMG_get_bc(params->per_mask(), &ipn);

		if (std::is_same<sten, seven_pt>::value) {
			BMG3_SymStd_SETUP_ITLI07_ex(fopd.data(), cop.data(), Pd.data(),
			                            fop.len(0), fop.len(1), fop.len(2),
			                            cop.len(0), cop.len(1), cop.len(2),
			                            ipn);
		} else {
			BMG3_SymStd_SETUP_ITLI27_ex(fopd.data(), cop.data(), Pd.data(),
			                            fop.len(0), fop.len(1), fop.len(2),
			                            cop.len(0), cop.len(1), cop.len(2),
			                            ipn);
		}
	}
};

}}

#endif
