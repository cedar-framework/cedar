#ifndef CEDAR_3D_MPI_PLANE_MPI_H
#define CEDAR_3D_MPI_PLANE_MPI_H

#include <vector>
#include <memory>
#include <map>

#include <mpi.h>
#ifdef PLANE_AGG
#include <abt.h>
#endif

#include <cedar/mpi/mpi_wrapper.h>

namespace cedar { namespace cdr3 { namespace mpi {

class plane_setup_mpi : public mpi_wrapper
{
public:
	plane_setup_mpi();
	plane_setup_mpi(std::vector<int> *keys);
	int comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) override;
	std::vector<int> *get_keys() { return key_data.get(); }

protected:
	bool ismaster;
	int currind;
	std::unique_ptr<std::vector<int>> key_data;
	std::vector<int> *keys;
	std::vector<std::unique_ptr<MPI_Comm>> comms;
};


#ifdef PLANE_AGG
class plane_mpi : public mpi_wrapper
{
public:
	plane_mpi(int nplanes, std::vector<ABT_thread> *threads);
	plane_mpi(int nplanes, int worker_id, std::vector<ABT_thread> *threads);
	int gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
	           void *recvbuf, int recvcount, MPI_Datatype recvtype,
	           int root, MPI_Comm comm) override;
	int gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
	            void *recvbuf, const int *recvcounts, const int *displs,
	            MPI_Datatype recvtype, int root, MPI_Comm comm);
	int scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
	            void *recvbuf, int recvcount, MPI_Datatype recvtype,
	            int root, MPI_Comm comm) override;
	int scatterv(const void *sendbuf, const int *sendcounts, const int *displs,
	             MPI_Datatype sendtype, void *recvbuf, int recvcount,
	             MPI_Datatype recvtype,
	             int root, MPI_Comm comm) override;

protected:
	int nplanes;
	bool ismaster;
	std::map<std::pair<int,int>, MPI_Datatype> tcache;
	std::vector<ABT_thread> *threads;
	int wid;

	MPI_Datatype get_aggtype(MPI_Comm comm, int plane_len, MPI_Datatype dtype);
};
#endif

}}}

#endif
