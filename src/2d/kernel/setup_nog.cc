#include <mpi.h>

#include "cedar/2d/kernel/setup_nog.h"

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_SETUP_nog(len_t NLx, len_t NLy, len_t NLXYc,
	                           int *NOG, len_t iGs, len_t jGs, int MPICOMM);
}


namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	using namespace cedar::cdr2;
	void fortran_setup_nog(const kernel_params & params, grid_topo & topo, len_t min_coarse, int *nog)
	{
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		BMG2_SymStd_SETUP_nog(topo.nlocal(0)-2, topo.nlocal(1)-2, min_coarse,
		                      nog, topo.is(0), topo.is(1), fcomm);
	}
}

}}}
