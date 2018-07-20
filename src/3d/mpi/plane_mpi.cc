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


}}}
