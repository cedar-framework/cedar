#include <cedar/2d/restrict.h>
#include <cedar/2d/ftn/BMG_parameters_c.h>

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_restrict(real_t*, real_t*, real_t*,
	                          int, int, int, int, int);
	void BMG_get_bc(int, int*);
}

using namespace cedar;
using namespace cedar::cdr2;

void restrict_f90::run(const restrict_op & R,
                       const grid_func & fine,
                       grid_func & coarse)
{
		int ibc;

		auto & fined = const_cast<grid_func&>(fine);
		auto & Rd = const_cast<restrict_op&>(R);
		prolong_op & P = Rd.getP();

		BMG_get_bc(params->per_mask(), &ibc);

		BMG2_SymStd_restrict(fined.data(), coarse.data(),
		                     P.data(), fined.len(0), fined.len(1),
		                     coarse.len(0), coarse.len(1), ibc);
}
