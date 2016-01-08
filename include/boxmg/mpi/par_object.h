#ifndef BOXMG_PAR_OBJECT_H
#define BOXMG_PAR_OBJECT_H

#include <mpi.h>
#include <boxmg/mpi/grid_topo.h>

namespace boxmg {

class par_object
{
public:
	MPI_Comm comm;
	grid_topo & grid() { return grid_; }
	const grid_topo & grid() const { return grid_; }

protected:
	topo_ptr grid_;
};

}


#endif
