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


	int gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
	           void *recvbuf, int recvcount, MPI_Datatype recvtype,
	           int root, MPI_Comm comm) override
	{
		return MPI_Gather(sendbuf, sendcount, sendtype,
		                  recvbuf, recvcount, recvtype,
		                  root, comm);
	}

	int gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
	                    void *recvbuf, const int *recvcounts, const int *displs,
	                    MPI_Datatype recvtype, int root, MPI_Comm comm) override
	{
		return MPI_Gatherv(sendbuf, sendcount, sendtype,
		                   recvbuf, recvcounts, displs,
		                   recvtype, root, comm);
	}

	int scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
	            void *recvbuf, int recvcount, MPI_Datatype recvtype,
	            int root, MPI_Comm comm) override
	{
		return MPI_Scatter(sendbuf, sendcount, sendtype,
		                   recvbuf, recvcount, recvtype,
		                   root, comm);
	}

	int scatterv(const void *sendbuf, const int *sendcounts, const int *displs,
	             MPI_Datatype sendtype, void *recvbuf, int recvcount,
	             MPI_Datatype recvtype,
	             int root, MPI_Comm comm) override
	{
		return MPI_Scatterv(sendbuf, sendcounts, displs,
		                    sendtype, recvbuf, recvcount,
		                    recvtype, root, comm);
	}
};

}

#endif
