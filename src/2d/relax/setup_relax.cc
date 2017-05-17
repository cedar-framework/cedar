#include "cedar/2d/ftn/BMG_parameters_c.h"
#include "cedar/2d/relax/setup_relax.h"

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_SETUP_recip(real_t *so, real_t *sor, len_t nx, len_t ny, int nstncl, int nsor_v);
	void BMG2_SymStd_SETUP_lines_x(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl, int JPN);
	void BMG2_SymStd_SETUP_lines_y(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl, int JPN);
	void BMG_get_bc(int, int*);
}


namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	using namespace cedar::cdr2;

	template<>
	void setup_rbgs_point(const kernel_params & params,
	                      const stencil_op<five_pt> & so,
	                      relax_stencil & sor)
	{
		len_t nx, ny;
		int nstencil, nsorv;

		auto & sod = const_cast<stencil_op<five_pt>&>(so);

		nx = so.len(0);
		ny = so.len(1);

		nstencil = 3;

		nsorv = 2;

		BMG2_SymStd_SETUP_recip(sod.data(), sor.data(), nx, ny, nstencil, nsorv);
	}

	template<>
	void setup_rbgs_point(const kernel_params & params,
	                      const stencil_op<nine_pt> & so,
	                      relax_stencil & sor)
	{
		len_t nx, ny;
		int nstencil, nsorv;

		auto & sod = const_cast<stencil_op<nine_pt>&>(so);

		nx = so.len(0);
		ny = so.len(1);

		nstencil = 5;

		nsorv = 2;

		BMG2_SymStd_SETUP_recip(sod.data(), sor.data(), nx, ny, nstencil, nsorv);
	}


	template<>
	void setup_rbgs_x(const kernel_params & params,
	                  const stencil_op<five_pt> & so,
	                  relax_stencil & sor)
	{
		int nx, ny, nstencil, jpn;

		auto & sod = const_cast<stencil_op<five_pt>&>(so);

		nx = so.len(0);
		ny = so.len(1);

		nstencil = 3;

		BMG_get_bc(params.per_mask(), &jpn);

		BMG2_SymStd_SETUP_lines_x(sod.data(), sor.data(), nx, ny, nstencil, jpn);
	}

	template<>
	void setup_rbgs_x(const kernel_params & params,
	                  const stencil_op<nine_pt> & so,
	                  relax_stencil & sor)
	{
		int nx, ny, nstencil, jpn;

		auto & sod = const_cast<stencil_op<nine_pt>&>(so);

		nx = so.len(0);
		ny = so.len(1);

		nstencil = 5;

		BMG_get_bc(params.per_mask(), &jpn);

		BMG2_SymStd_SETUP_lines_x(sod.data(), sor.data(), nx, ny, nstencil, jpn);
	}


	template<>
	void setup_rbgs_y(const kernel_params & params,
	                  const stencil_op<five_pt> & so,
	                  relax_stencil & sor)
	{
		int nx, ny, nstencil, jpn;

		auto & sod = const_cast<stencil_op<five_pt>&>(so);

		nx = so.len(0);
		ny = so.len(1);

		nstencil = 3;

		BMG_get_bc(params.per_mask(), &jpn);

		BMG2_SymStd_SETUP_lines_y(sod.data(), sor.data(), nx, ny, nstencil, jpn);
	}


	template<>
	void setup_rbgs_y(const kernel_params & params,
	                  const stencil_op<nine_pt> & so,
	                  relax_stencil & sor)
	{
		int nx, ny, nstencil, jpn;

		auto & sod = const_cast<stencil_op<nine_pt>&>(so);

		nx = so.len(0);
		ny = so.len(1);

		nstencil = 5;

		BMG_get_bc(params.per_mask(), &jpn);

		BMG2_SymStd_SETUP_lines_y(sod.data(), sor.data(), nx, ny, nstencil, jpn);
	}
}

}}}
