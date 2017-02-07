#ifndef BOXMG_PAR_OBJECT_H
#define BOXMG_PAR_OBJECT_H

#include <mpi.h>
#include <boxmg/mpi/grid_topo.h>

namespace boxmg {

class par_object
{
public:
	par_object() {}
par_object(topo_ptr grd) : grid_(grd) {}
par_object(topo_ptr grd, MPI_Comm comm) : comm(comm), grid_(grd) {}
	MPI_Comm comm;
	grid_topo & grid() { return *grid_; }
	const grid_topo & grid() const { return *grid_; }
	topo_ptr grid_ptr() const { return grid_; }
	void *halo_ctx;

protected:
	topo_ptr grid_;
};

}


#endif
