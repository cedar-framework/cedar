#include "boxmg/2d/ftn/mpi/BMG_workspace_c.h"
#include "boxmg/2d/ftn/BMG_parameters_c.h"

#include <boxmg/3d/cg/setup_cg_lu.h>

extern "C" {
	using namespace boxmg;
	void BMG3_SymStd_SETUP_cg_LU(real_t *so, len_t ii, len_t jj, len_t kk, int NStncl,
	                             real_t *abd, len_t nabd1, len_t nabd2, int ibc);

}


namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void setup_cg_lu(const stencil_op &so,
	                 grid_func &ABD)
	{
		int nstencil, ibc;
		auto & sod = const_cast<stencil_op&>(so);
		auto & so_sten = sod.stencil();

		ibc = BMG_BCs_definite;

		if (so_sten.five_pt()) {
			nstencil = 4;
		} else {
			nstencil = 14;
		}

		BMG3_SymStd_SETUP_cg_LU(so_sten.data(), so_sten.len(0), so_sten.len(1), so_sten.len(2),
		                        nstencil, ABD.data(), ABD.len(0), ABD.len(1), ibc);
	}

}

}}}
