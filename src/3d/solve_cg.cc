#include <cedar/3d/solve_cg.h>


extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SOLVE_cg(real_t *q, real_t *qf,
	                          len_t ii, len_t jj, len_t kk,
	                          real_t *abd, real_t *bbd, len_t nabd1, len_t nabd2,
	                          int ibc);
	void BMG_get_bc(int, int*);
}


namespace cedar { namespace cdr3 {

void solve_cg_f90::run(grid_func & x,
                       const grid_func & b,
                       const grid_func & ABD,
                       real_t *bbd)
{
		int ibc;

		auto & bd = const_cast<grid_func&>(b);
		auto & abd_data = const_cast<grid_func&>(ABD);

		BMG_get_bc(params->per_mask(), &ibc);

		BMG3_SymStd_SOLVE_cg(x.data(), bd.data(), x.len(0), x.len(1), x.len(2),
		                     abd_data.data(), &bbd[0], ABD.len(0), ABD.len(1), ibc);
}

}}
