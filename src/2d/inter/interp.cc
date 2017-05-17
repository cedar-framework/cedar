#include "cedar/2d/ftn/BMG_parameters_c.h"
#include "cedar/2d/inter/interp.h"

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_interp_add(real_t*, real_t*, real_t*,real_t*,real_t*,
	                            len_t, len_t, len_t, len_t, int, int);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	void fortran_interp(const kernel_params & params,
	                    const inter::prolong_op & P,
	                    const grid_func & coarse,
	                    const grid_func & residual,
	                    grid_func & fine)
	{
		using namespace cedar::cdr2;

		int nstencil, ibc;

		inter::prolong_op & Pd = const_cast<inter::prolong_op&>(P);
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

		BMG_get_bc(params.per_mask(), &ibc);

		BMG2_SymStd_interp_add(fine.data(), coarsed.data(), res.data(), fop_data, Pd.data(),
		                       coarsed.len(0), coarsed.len(1), fine.len(0), fine.len(1),
		                       nstencil, ibc);
	}
}


}}}
