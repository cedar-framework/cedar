#include "cedar/2d/cg/setup_cg_lu.h"

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_SETUP_cg_LU(real_t*, len_t*, len_t*, int*, real_t*, len_t*,len_t*,int*);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	using namespace cedar::cdr2;

	template<>
	void setup_cg_lu(const kernel_params & params, const stencil_op<five_pt> & so,
	                 grid_func & ABD)
	{
		len_t nx, ny;
		int nstencil;
		len_t nabd1, nabd2;
		int ibc;

		auto & sod = const_cast<stencil_op<five_pt>&>(so);

		nx = so.len(0);
		ny = so.len(1);

		nstencil = 3;

		nabd1 = ABD.len(0);
		nabd2 = ABD.len(1);

		BMG_get_bc(params.per_mask(), &ibc);

		BMG2_SymStd_SETUP_cg_LU(sod.data(), &nx, &ny, &nstencil,
		                        ABD.data(), &nabd1, &nabd2, &ibc);
	}

	template<>
	void setup_cg_lu(const kernel_params & params, const stencil_op<nine_pt> & so,
	                 grid_func & ABD)
	{
		len_t nx, ny;
		int nstencil;
		len_t nabd1, nabd2;
		int ibc;

		auto & sod = const_cast<stencil_op<nine_pt>&>(so);

		nx = so.len(0);
		ny = so.len(1);

		nstencil = 5;

		nabd1 = ABD.len(0);
		nabd2 = ABD.len(1);

		BMG_get_bc(params.per_mask(), &ibc);

		BMG2_SymStd_SETUP_cg_LU(sod.data(), &nx, &ny, &nstencil,
		                        ABD.data(), &nabd1, &nabd2, &ibc);
	}
}

}}}
