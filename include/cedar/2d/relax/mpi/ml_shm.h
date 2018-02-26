#ifndef CEDAR_2D_KERNEL_RELAX_MPI_ML_SHM_H
#define CEDAR_2D_KERNEL_RELAX_MPI_ML_SHM_H

#include <mpi.h>
#include <cedar/types.h>

namespace cedar { namespace cdr2 { namespace mpi {

class ml_relax_pup
{
public:
	void pack(real_t *rwork, int nlines);
	void unpack(real_t *rwork, int nlines);
	void init(MPI_Comm shm_comm, MPI_Comm node_comm,
	          MPI_Win shm_win, real_t *shm_buff);

protected:
	MPI_Comm shm_comm;
	MPI_Comm node_comm;
	MPI_Win shm_win;
	real_t *shm_buff;
};

}}}

#endif
