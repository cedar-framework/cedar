#ifndef CEDAR_3D_SETUP_RELAX_MPI_H
#define CEDAR_3D_SETUP_RELAX_MPI_H

#include <cedar/kernel_params.h>
#include <cedar/2d/ftn/BMG_parameters_c.h>
#include <cedar/3d/mpi/stencil_op.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG3_SymStd_SETUP_recip(real_t *so, real_t *sor,
	                                 len_t nx, len_t ny, len_t nz,
	                                 int nstencil, int nsorv);
}

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;

	template<class sten>
	void mpi_setup_rbgs_point(const kernel_params & params,
	                          const mpi::stencil_op<sten> & so,
	                          relax_stencil & sor)
	{
		int nstencil, nsorv;

		auto & sod = const_cast<mpi::stencil_op<sten>&>(so);

		nstencil = stencil_ndirs<sten>::value;

		nsorv = 2;

		MPI_BMG3_SymStd_SETUP_recip(sod.data(),
		                            sor.data(),
		                            so.len(0), so.len(1), so.len(2),
		                            nstencil, nsorv);
	}
}

}}}


#endif
