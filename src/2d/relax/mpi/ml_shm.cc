#include <cstring>

#include <cedar/util/log.h>
#include <cedar/2d/relax/mpi/ml_shm.h>

namespace cedar { namespace cdr2 { namespace mpi {

void ml_relax_pup::pack(real_t *rwork, int nlines)
{
	MPI_Win_lock_all(MPI_MODE_NOCHECK, shm_win);

	std::memcpy(shm_buff, rwork, nlines*8*sizeof(real_t));

	MPI_Win_sync(shm_win);
	MPI_Barrier(shm_comm);

	int rank, size;
	MPI_Comm_rank(shm_comm, &rank);
	MPI_Comm_size(shm_comm, &size);

	if (rank == 0) {
		real_t *seg;
		MPI_Aint seg_size;
		int disp_unit;
		for (int proc = 1; proc < size; proc++) {
			MPI_Win_shared_query(shm_win, proc, &seg_size, &disp_unit, &seg);
			std::memcpy(rwork + (proc * nlines * 8), seg, nlines*8*sizeof(real_t));
		}
	}

	MPI_Win_unlock_all(shm_win);
}


void ml_relax_pup::unpack(real_t *rwork, int nlines)
{
	int rank, size;
	MPI_Comm_rank(shm_comm, &rank);
	MPI_Comm_size(shm_comm, &size);

	MPI_Win_lock_all(MPI_MODE_NOCHECK, shm_win);
	if (rank == 0) {
		real_t *seg;
		MPI_Aint seg_size;
		int disp_unit;
		for (int proc = 1; proc < size; proc++) {
			MPI_Win_shared_query(shm_win, proc, &seg_size, &disp_unit, &seg);
			std::memcpy(seg, rwork + (proc * nlines * 8), nlines*8*sizeof(real_t));
		}
	}

	MPI_Win_sync(shm_win);
	MPI_Barrier(shm_comm);

	if (rank != 0)
		std::memcpy(rwork, shm_buff, nlines * 8 * sizeof(real_t));

	MPI_Win_unlock_all(shm_win);
}


void ml_relax_pup::init(MPI_Comm shm_comm, MPI_Comm node_comm,
                        MPI_Win shm_win, real_t *shm_buff)
{
	this->shm_comm = shm_comm;
	this->node_comm = node_comm;
	this->shm_win = shm_win;
	this->shm_buff = shm_buff;
}

}}}
