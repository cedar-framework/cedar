#ifndef CEDAR_3D_MPI_PLANE_MPI_H
#define CEDAR_3D_MPI_PLANE_MPI_H

#include <vector>
#include <memory>
#include <map>

#include <mpi.h>
#include <abt.h>

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


class plane_mpi : public mpi_wrapper
{
public:
	plane_mpi(int nplanes);
	plane_mpi(int nplanes, ABT_barrier barrier);
	int gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
	           void *recvbuf, int recvcount, MPI_Datatype recvtype,
	           int root, MPI_Comm comm) override;
	int scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
	            void *recvbuf, int recvcount, MPI_Datatype recvtype,
	            int root, MPI_Comm comm) override;

	ABT_barrier get_barrier() { return barrier; }

protected:
	int nplanes;
	bool ismaster;
	std::map<std::pair<int,int>, MPI_Datatype> tcache;
	ABT_barrier barrier;

	MPI_Datatype get_aggtype(MPI_Comm comm, int plane_len, MPI_Datatype dtype);
};

}}}

#endif
