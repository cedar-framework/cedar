#include <cedar/util/log.h>
#include <cedar/mpi/mpi_wrapper.h>

namespace cedar {

class master_splitter : public mpi_wrapper
{
public:
	master_splitter(std::vector<int> *keys): keys(keys) {}
	int comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) override
	{
		int kvkey;
		MPI_Comm_create_keyval(MPI_NULL_COPY_FN,
		                       MPI_NULL_DELETE_FN,
		                       &kvkey, (void*) 0);
		keys->push_back(kvkey);
		MPI_Comm_split(comm, color, key, newcomm);
		MPI_Comm_set_attr(comm, kvkey, newcomm);

		return 0;
	}

protected:
	std::vector<int> *keys;
	std::vector<int> testvals;
};


class worker_splitter : public mpi_wrapper
{
public:
	worker_splitter(std::vector<int> *keys) : currind(0), keys(keys) {}
	int comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) override
	{
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
		return 0;
	}

protected:
	int currind;
	std::vector<int> *keys;
};

}
