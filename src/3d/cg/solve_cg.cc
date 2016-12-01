#include "boxmg/2d/ftn/BMG_parameters_c.h"

#include <boxmg/3d/cg/solve_cg.h>


extern "C" {
	using namespace boxmg;
	void BMG3_SymStd_SOLVE_cg(real_t *q, real_t *qf,
	                          len_t ii, len_t jj, len_t kk,
	                          real_t *abd, real_t *bbd, len_t nabd1, len_t nabd2,
	                          int ibc);
}

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void fortran_solve_cg(grid_func & x,
	                      const grid_func & b,
	                      const grid_func & ABD,
	                      real_t *bbd)
	{
		using namespace boxmg::bmg3;
		int ibc;

		auto & bd = const_cast<grid_func&>(b);
		auto & abd_data = const_cast<grid_func&>(ABD);

		ibc = BMG_BCs_definite;

		BMG3_SymStd_SOLVE_cg(x.data(), bd.data(), x.len(0), x.len(1), x.len(2),
		                     abd_data.data(), &bbd[0], ABD.len(0), ABD.len(1), ibc);
	}
}

}}}

