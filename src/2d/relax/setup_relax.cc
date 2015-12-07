#include "boxmg/2d/ftn/BMG_parameters_c.h"
#include "boxmg/2d/relax/setup_relax.h"

extern "C" {
	void bmg2_symstd_setup_recip(double*, double*, int*, int*, int*, int*);
	using namespace boxmg;
	void BMG2_SymStd_SETUP_lines_x(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl, int JPN);
	void BMG2_SymStd_SETUP_lines_y(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl, int JPN);
	void MPI_BMG2_SymStd_SETUP_recip(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
	void MPI_BMG2_SymStd_SETUP_lines_x(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
	void MPI_BMG2_SymStd_SETUP_lines_y(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
}


namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	using namespace boxmg::bmg2d;
	void setup_rbgs_point(const stencil_op & so,
	                      relax_stencil & sor)
	{
		int nx, ny, nstencil, nsorv;

		const grid_stencil & so_sten = so.stencil();
		stencil_op & sod = const_cast<stencil_op&>(so);

		nx = so_sten.len(0);
		ny = so_sten.len(1);

		if (so_sten.five_pt()) nstencil = 3;
		else nstencil = 5;

		nsorv = 2;

		bmg2_symstd_setup_recip(sod.data(), sor.data(), &nx, &ny, &nstencil, &nsorv);
	}


	void setup_rbgs_x(const stencil_op & so,
	                  relax_stencil & sor)
	{
		int nx, ny, nstencil, jpn;

		const grid_stencil & so_sten = so.stencil();
		stencil_op & sod = const_cast<stencil_op&>(so);

		nx = so_sten.len(0);
		ny = so_sten.len(1);

		if (so_sten.five_pt()) nstencil = 3;
		else nstencil = 5;

		jpn = BMG_BCs_definite;

		BMG2_SymStd_SETUP_lines_x(sod.data(), sor.data(), nx, ny, nstencil, jpn);
	}


	void setup_rbgs_y(const stencil_op & so,
	                  relax_stencil & sor)
	{
		int nx, ny, nstencil, jpn;

		const grid_stencil & so_sten = so.stencil();
		stencil_op & sod = const_cast<stencil_op&>(so);

		nx = so_sten.len(0);
		ny = so_sten.len(1);

		if (so_sten.five_pt()) nstencil = 3;
		else nstencil = 5;

		jpn = BMG_BCs_definite;

		BMG2_SymStd_SETUP_lines_y(sod.data(), sor.data(), nx, ny, nstencil, jpn);
	}



	void mpi_setup_rbgs_point(const stencil_op & so,
	                          relax_stencil & sor)
	{
		int nx, ny, nstencil;

		const grid_stencil & so_sten = so.stencil();
		stencil_op & sod = const_cast<stencil_op&>(so);

		nx = so_sten.len(0);
		ny = so_sten.len(1);

		if (so_sten.five_pt()) nstencil = 3;
		else nstencil = 5;

		MPI_BMG2_SymStd_SETUP_recip(sod.data(), sor.data(), nx, ny, nstencil);
	}


	void mpi_setup_rbgs_x(const stencil_op & so,
	                      relax_stencil & sor)
	{
		int nx, ny, nstencil;

		const grid_stencil & so_sten = so.stencil();
		stencil_op & sod = const_cast<stencil_op&>(so);

		nx = so_sten.len(0);
		ny = so_sten.len(1);

		if (so_sten.five_pt()) nstencil = 3;
		else nstencil = 5;

		MPI_BMG2_SymStd_SETUP_lines_x(sod.data(), sor.data(), nx, ny, nstencil);
	}


	void mpi_setup_rbgs_y(const stencil_op & so,
	                      relax_stencil & sor)
	{
		int nx, ny, nstencil;

		const grid_stencil & so_sten = so.stencil();
		stencil_op & sod = const_cast<stencil_op&>(so);

		nx = so_sten.len(0);
		ny = so_sten.len(1);

		if (so_sten.five_pt()) nstencil = 3;
		else nstencil = 5;

		MPI_BMG2_SymStd_SETUP_lines_y(sod.data(), sor.data(), nx, ny, nstencil);
	}
}

}}}
