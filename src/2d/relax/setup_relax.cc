#include "cedar/2d/ftn/BMG_parameters_c.h"
#include "cedar/2d/relax/setup_relax.h"

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_SETUP_recip(real_t *so, real_t *sor, len_t nx, len_t ny, int nstncl, int nsor_v);
	void BMG2_SymStd_SETUP_lines_x(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl, int JPN);
	void BMG2_SymStd_SETUP_lines_y(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl, int JPN);
}


namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	using namespace cedar::cdr2;
	void setup_rbgs_point(const kernel_params & params,
	                      const stencil_op & so,
	                      relax_stencil & sor)
	{
		len_t nx, ny;
		int nstencil, nsorv;

		const grid_stencil & so_sten = so.stencil();
		stencil_op & sod = const_cast<stencil_op&>(so);

		nx = so_sten.len(0);
		ny = so_sten.len(1);

		if (so_sten.five_pt()) nstencil = 3;
		else nstencil = 5;

		nsorv = 2;

		BMG2_SymStd_SETUP_recip(sod.data(), sor.data(), nx, ny, nstencil, nsorv);
	}


	void setup_rbgs_x(const kernel_params & params,
	                  const stencil_op & so,
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


	void setup_rbgs_y(const kernel_params & params,
	                  const stencil_op & so,
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
}

}}}
