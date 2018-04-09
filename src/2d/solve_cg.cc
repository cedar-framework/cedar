#include <cedar/2d/solve_cg.h>
#include <cedar/2d/ftn/BMG_parameters_c.h>

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_SOLVE_cg(real_t*, real_t*, len_t, len_t, real_t*, real_t*, len_t, len_t, int);
	void BMG_get_bc(int, int*);
}


namespace cedar { namespace cdr2 {

void solve_cg_f90::run(grid_func & x,
                       const grid_func & b,
                       const grid_func & ABD,
                       real_t * bbd)
{
		int ibc;

		grid_func & bd = const_cast<grid_func&>(b);
		grid_func & abd_data = const_cast<grid_func&>(ABD);


		BMG_get_bc(params->per_mask(), &ibc);

		BMG2_SymStd_SOLVE_cg(x.data(), bd.data(), x.len(0), x.len(1),
		                    abd_data.data(), &bbd[0], ABD.len(0), ABD.len(1), ibc);
}
}}
