#include <boxmg/2d/ftn/BMG_parameters_c.h>

#include <boxmg/3d/inter/restrict.h>

extern "C" {
	using namespace boxmg;
	void BMG3_SymStd_restrict(real_t *q, real_t *qc, real_t *ci,
	                          len_t nx, len_t ny, len_t nz,
	                          len_t nxc, len_t nyc, len_t nzc,
	                          int jpn);
}

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void fortran_restrict(const inter::restrict_op & R,
	                      const grid_func & fine,
	                      grid_func & coarse)
	{
		using namespace boxmg::bmg3;
		int ibc;

		auto & fined = const_cast<grid_func&>(fine);
		auto & Rd = const_cast<inter::restrict_op&>(R);
		inter::prolong_op &P = Rd.getP();
		ibc = BMG_BCs_definite;

		BMG3_SymStd_restrict(fined.data(), coarse.data(),
		                     P.data(),
		                     fined.len(0), fined.len(1), fined.len(2),
		                     coarse.len(0), coarse.len(1), coarse.len(2),
		                     ibc);
	}
}

}}}
