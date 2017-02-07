#include "boxmg/2d/ftn/BMG_parameters_c.h"

#include "boxmg/2d/inter/restrict.h"

extern "C" {
	using namespace boxmg;
	void BMG2_SymStd_restrict(real_t*, real_t*, real_t*,
	                          int, int, int, int, int);
}

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void fortran_restrict(const inter::restrict_op & R,
	                      const grid_func & fine,
	                      grid_func & coarse)
	{
		using namespace boxmg::bmg2d;
		int ibc;

		auto & fined = const_cast<grid_func&>(fine);
		auto & Rd = const_cast<inter::restrict_op&>(R);
		inter::prolong_op & P = Rd.getP();
		ibc = BMG_BCs_definite;

		BMG2_SymStd_restrict(fined.data(), coarse.data(),
		                     P.data(), fined.len(0), fined.len(1),
		                     coarse.len(0), coarse.len(1), ibc);
	}
}

}}}
