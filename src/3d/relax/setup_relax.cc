#include <cedar/2d/ftn/BMG_parameters_c.h>
#include <cedar/3d/relax/setup_relax.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SETUP_recip(real_t *so, real_t *sor,
	                             len_t nx, len_t ny, len_t nz,
	                             int nstencl, int nsorv);
}


namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	using namespace cedar::cdr3;
	void setup_rbgs_point(const kernel_params & params,
	                      const stencil_op &so,
	                      relax_stencil &sor)
	{
		int nsorv, nstencil;
		const grid_stencil &sten = so.stencil();
		stencil_op &sod = const_cast<stencil_op&>(so);

		if (sten.five_pt()) {
			nstencil = 4;
		} else {
			nstencil = 14;
		}

		nsorv = 2;

		BMG3_SymStd_SETUP_recip(sod.data(),
		                        sor.data(),
		                        sten.len(0), sten.len(1), sten.len(2),
		                        nstencil, nsorv);
	}
}

}}}
