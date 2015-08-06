#include <mpi.h>

#include "setup_nog.h"

extern "C" {
	using namespace boxmg;
	void BMG2_SymStd_SETUP_nog(len_t NLx, len_t NLy, len_t NLXYc,
	                           int *NOG, len_t iGs, len_t jGs, int MPICOMM);
}


namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	using namespace boxmg::bmg2d::core;
	void fortran_setup_nog(mpi::GridTopo & topo, len_t min_coarse, int *nog)
	{
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		BMG2_SymStd_SETUP_nog(topo.nlocal(0), topo.nlocal(1), min_coarse,
		                      nog, topo.is(0), topo.is(1), fcomm);
	}
}

}}}
