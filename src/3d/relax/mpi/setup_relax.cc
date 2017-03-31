#include <cedar/2d/ftn/BMG_parameters_c.h>
#include <cedar/3d/relax/setup_relax.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG3_SymStd_SETUP_recip(real_t *so, real_t *sor,
	                                 len_t nx, len_t ny, len_t nz,
	                                 int nstencil, int nsorv);
}


namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	using namespace cedar::cdr3;
	void mpi_setup_rbgs_point(const mpi::stencil_op & so,
	                          relax_stencil & sor)
	{
		int nstencil, nsorv;

		auto & so_sten = so.stencil();
		auto & sod = const_cast<mpi::stencil_op&>(so);

		if (so_sten.five_pt()) nstencil = 4;
		else nstencil = 14;

		nsorv = 2;

		MPI_BMG3_SymStd_SETUP_recip(sod.data(),
		                            sor.data(),
		                            so_sten.len(0), so_sten.len(1), so_sten.len(2),
		                            nstencil, nsorv);
	}
}

}}}
