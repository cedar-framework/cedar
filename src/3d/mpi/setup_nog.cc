#include <cedar/3d/mpi/setup_nog.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SETUP_nog(len_t nlx, len_t nly, len_t nlz, len_t nlxyzc,
	                           int *nog, len_t igs, len_t jgs, len_t kgs, int MPICOMM);
}

using namespace cedar;
using namespace cedar::cdr3;

void mpi::setup_nog_f90::run(grid_topo & topo, len_t min_coarse, int *nog)
{
	MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

	BMG3_SymStd_SETUP_nog(topo.nlocal(0)-2, topo.nlocal(1)-2, topo.nlocal(2)-2, min_coarse,
	                      nog, topo.is(0), topo.is(1), topo.is(2), fcomm);
}
