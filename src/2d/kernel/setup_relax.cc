#include "setup_relax.h"

extern "C" {
	void bmg2_symstd_setup_recip(double*, double*, int*, int*, int*, int*);
	using namespace boxmg;
	void MPI_BMG2_SymStd_SETUP_recip(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
}


namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	using namespace boxmg::bmg2d::core;
	void setup_rbgs_point(const StencilOp & so,
	                      RelaxStencil & sor)
	{
		int nx, ny, nstencil, nsorv;

		const GridStencil & so_sten = so.stencil();
		core::StencilOp & sod = const_cast<core::StencilOp&>(so);

		nx = so_sten.len(0);
		ny = so_sten.len(1);

		if (so_sten.five_pt()) nstencil = 3;
		else nstencil = 5;

		nsorv = 2;

		bmg2_symstd_setup_recip(sod.data(), sor.data(), &nx, &ny, &nstencil, &nsorv);
	}


	void mpi_setup_rbgs_point(const StencilOp & so,
	                          RelaxStencil & sor)
	{
		int nx, ny, nstencil;

		const GridStencil & so_sten = so.stencil();
		core::StencilOp & sod = const_cast<core::StencilOp&>(so);

		nx = so_sten.len(0);
		ny = so_sten.len(1);

		if (so_sten.five_pt()) nstencil = 3;
		else nstencil = 5;

		MPI_BMG2_SymStd_SETUP_recip(sod.data(), sor.data(), nx, ny, nstencil);
	}
}

}}}
