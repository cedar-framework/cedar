#include <cassert>
#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"

#include <cedar/array.h>
#include <cedar/mpi/block_partition.h>
#include <cedar/2d/mpi/redist_solver.h>


extern "C" {
	void MSG_pause(MPI_Fint *msg_comm);
	void MSG_play(MPI_Fint msg_comm);
}

using namespace cedar;
using namespace cedar::cdr2::mpi;





