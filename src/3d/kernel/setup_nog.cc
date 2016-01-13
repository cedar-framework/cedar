#include <mpi.h>

#include "boxmg/3d/kernel/setup_nog.h"

extern "C" {
	using namespace boxmg;
	void BMG3_SymStd_SETUP_nog(len_t nlx, len_t nly, len_t nlz, len_t nlxyzc,
	                           int *nog, len_t igs, len_t jgs, len_t kgs, int MPICOMM);
}


namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	using namespace boxmg::bmg3;
	void fortran_setup_nog(grid_topo & topo, len_t min_coarse, int *nog)
	{
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		BMG3_SymStd_SETUP_nog(topo.nlocal(0), topo.nlocal(1), topo.nlocal(2), min_coarse,
		                      nog, topo.is(0), topo.is(1), topo.is(2), fcomm);
	}
}

}}}
