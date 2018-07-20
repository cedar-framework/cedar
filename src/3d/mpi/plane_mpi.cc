#include <cedar/types.h>
#include <cedar/3d/mpi/plane_mpi.h>

namespace cedar { namespace cdr3 { namespace mpi {

plane_setup_mpi::plane_setup_mpi() : ismaster(true), currind(0)
{
	key_data = std::make_unique<std::vector<int>>();
	keys = key_data.get();
}

plane_setup_mpi::plane_setup_mpi(std::vector<int> * keys) :
	ismaster(false), currind(0), keys(keys) {}

int plane_setup_mpi::comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
{
	if (ismaster) {
		int kvkey;
		MPI_Comm_create_keyval(MPI_NULL_COPY_FN,
		                       MPI_NULL_DELETE_FN,
		                       &kvkey, (void*) 0);
		keys->push_back(kvkey);
		MPI_Comm_split(comm, color, key, newcomm);
		comms.push_back(*newcomm);
		MPI_Comm_set_attr(comm, kvkey, &comms.back());
	} else {
		if (currind >= keys->size()) {
			log::error << "master has not split communicator" << std::endl;
			return 1;
		}

		int flag;
		MPI_Comm *master_comm;
		MPI_Comm_get_attr(comm, (*keys)[currind], &master_comm, &flag);
		*newcomm = *master_comm;
		currind++;

		if (not flag)
			log::error << "com_split attribute not found" << std::endl;
	}

	return 0;
}


MPI_Datatype plane_mpi::get_aggtype(MPI_Comm comm, int plane_len, MPI_Datatype dtype)
{
	int nproc;
	MPI_Comm_size(comm, &nproc);
	int type_size;
	MPI_Type_size(dtype, &type_size);

	auto atypeit = tcache.find(std::make_pair(nproc, plane_len));
	if (atypeit != tcache.end())
		return atypeit->second;

	MPI_Datatype plane_type, agg_type;
	int *displs = new int[nplanes];
	for (std::size_t i = 0; i < nplanes; i++)
		displs[i] = i*nproc*plane_len;
	MPI_Type_create_indexed_block(nplanes, plane_len, displs, dtype, &plane_type);
	MPI_Type_commit(&plane_type);
	MPI_Type_create_resized(plane_type, 0, plane_len*type_size, &agg_type);
	MPI_Type_commit(&agg_type);

	delete[] displs;

	tcache[std::make_pair(nproc, plane_len)] = agg_type;
	return agg_type;
}


int plane_mpi::gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                      void *recvbuf, int recvcount, MPI_Datatype recvtype,
                      int root, MPI_Comm comm)
{
	if (not ismaster)
		return 0;

	auto agg_type = get_aggtype(comm, sendcount, sendtype);

	int ierr = MPI_Gather(sendbuf, sendcount*nplanes, sendtype,
	                      recvbuf, 1, agg_type, 0, comm);

	return ierr;
}


int plane_mpi::scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                       void *recvbuf, int recvcount, MPI_Datatype recvtype,
                       int root, MPI_Comm comm)
{
	if (not ismaster)
		return 0;

	auto agg_type = get_aggtype(comm, sendcount, sendtype);

	int ierr = MPI_Scatter(sendbuf, 1, agg_type,
	                       recvbuf, sendcount*nplanes, recvtype, 0, comm);

	return ierr;
}


}}}
