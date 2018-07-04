#ifndef CEDAR_MPI_WRAPPER_H
#define CEDAR_MPI_WRAPPER_H

#include <mpi.h>

#include <cedar/services/message_passing.h>

namespace cedar {

class mpi_wrapper : public services::message_passing
{
public:
	int comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) override
	{
		return MPI_Comm_split(comm, color, key, newcomm);
	}
};

}

#endif
