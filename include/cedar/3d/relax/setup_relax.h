#ifndef CEDAR_3D_SETUP_RELAX_H
#define CEDAR_3D_SETUP_RELAX_H

#include <cedar/kernel_params.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/relax_stencil.h>
//#include <cedar/3d/relax/setup_planes.h>
#include <cedar/2d/ftn/BMG_parameters_c.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SETUP_recip(real_t *so, real_t *sor,
	                             len_t nx, len_t ny, len_t nz,
	                             int nstencl, int nsorv);
}

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	template<class sten>
	void setup_rbgs_point(const kernel_params & params,
	                      const stencil_op<sten> & so,
	                      relax_stencil & sor)
	{
		int nsorv, nstencil;
		auto &sod = const_cast<stencil_op<sten>&>(so);

		nstencil = stencil_ndirs<sten>::value;

		nsorv = 2;

		BMG3_SymStd_SETUP_recip(sod.data(),
		                        sor.data(),
		                        so.len(0), so.len(1), so.len(2),
		                        nstencil, nsorv);
	}
}

}}}

#endif
