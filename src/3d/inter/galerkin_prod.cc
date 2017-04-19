#include <cedar/2d/ftn/BMG_parameters_c.h>
#include <cedar/3d/inter/galerkin_prod.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SETUP_ITLI07_ex(real_t *so, real_t *soc, real_t *ci, len_t iif, len_t jjf,
	                                 len_t kkf, len_t iic, len_t jjc, len_t kkc,
	                                 int ipn);
	void BMG3_SymStd_SETUP_ITLI27_ex(real_t *so, real_t *soc, real_t *ci, len_t iif, len_t jjf,
	                                 len_t kkf, len_t iic, len_t jjc, len_t kkc, int ipn);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	using namespace cedar::cdr3;

	void galerkin_prod(const kernel_params & params,
	                   int kf, int kc, int nog,
	                   const inter::prolong_op &P,
	                   const stencil_op & fop,
	                   stencil_op & cop)
	{
		int ipn;
		const grid_stencil &fsten = fop.stencil();
		grid_stencil &csten = cop.stencil();
		stencil_op &fopd = const_cast<stencil_op&>(fop);
		inter::prolong_op &Pd = const_cast<inter::prolong_op&>(P);

		BMG_get_bc(params.per_mask(), &ipn);

		if (fsten.five_pt()) {
			BMG3_SymStd_SETUP_ITLI07_ex(fopd.data(), cop.data(), Pd.data(),
			                            fsten.len(0), fsten.len(1), fsten.len(2),
			                            csten.len(0), csten.len(1), csten.len(2),
			                            ipn);
		} else {
			BMG3_SymStd_SETUP_ITLI27_ex(fopd.data(), cop.data(), Pd.data(),
			                            fsten.len(0), fsten.len(1), fsten.len(2),
			                            csten.len(0), csten.len(1), csten.len(2),
			                            ipn);
		}
	}
}

}}}
