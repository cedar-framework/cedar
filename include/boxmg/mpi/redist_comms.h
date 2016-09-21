#ifndef BOXMG_2D_MPI_REDIST_COMMS_H
#define BOXMG_2D_MPI_REDIST_COMMS_H

#include <mpi.h>

struct redist_comms
{
	// Collective communicator for gathering data in processor blocks
	MPI_Comm pblock_comm;
	 /** Collective communicator for scattering data in active processors within a
	     processor block */
	MPI_Comm active_pblock_comm;
	MPI_Comm parent_comm; /** Parent solver's communicator */
	MPI_Comm redist_comm; /** Redistributed solver's communicator */
};

#endif
