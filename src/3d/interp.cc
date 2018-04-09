#include <cedar/2d/ftn/BMG_parameters_c.h>

#include <cedar/3d/interp.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_interp_add(real_t *q, real_t *qc,
	                            real_t *so, real_t *res, real_t *ci,
	                            len_t iic, len_t jjc, len_t kkc,
	                            len_t iif, len_t jjf, len_t kkf,
	                            int NStncl, int jpn);
	void BMG_get_bc(int, int*);
}


namespace cedar { namespace cdr3 {

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

		if (Pd.fine_is_seven) {
			nstencil = 4;
			fop_data = Pd.fine_op_seven->data();
		} else {
			nstencil = 14;
			fop_data = Pd.fine_op_xxvii->data();
		}

		BMG_get_bc(params->per_mask(), &ibc);

		BMG3_SymStd_interp_add(fine.data(), coarsed.data(), fop_data, res.data(), Pd.data(),
		                       coarsed.len(0), coarsed.len(1), coarsed.len(2),
		                       fine.len(0), fine.len(1), fine.len(2),
		                       nstencil, ibc);
}

}}
