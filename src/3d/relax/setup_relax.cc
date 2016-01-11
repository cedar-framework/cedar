#include <boxmg/2d/ftn/BMG_parameters_c.h>
#include <boxmg/3d/relax/setup_relax.h>

extern "C" {
	using namespace boxmg;
	void BMG3_SymStd_SETUP_recip(real_t *so, real_t *sor,
	                             len_t nx, len_t ny, len_t nz,
	                             int nstencl, int nsorv);
	void MPI_BMG3_SymStd_SETUP_recip(real_t *so, real_t *sor,
	                                 len_t nx, len_t ny, len_t nz,
	                                 int nstencil, int nsorv);
}


namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	using namespace boxmg::bmg3;
	void setup_rbgs_point(const stencil_op &so,
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
