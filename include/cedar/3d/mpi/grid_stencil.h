#ifndef CEDAR_3D_MPI_GRID_STENCIL_H
#define CEDAR_3D_MPI_GRID_STENCIL_H

#include <mpi.h>
#include <cedar/3d/grid_stencil.h>
#include <cedar/mpi/grid_topo.h>

namespace cedar { namespace cdr3 { namespace mpi {

class grid_stencil : public cdr3::grid_stencil
{
public:
	grid_stencil();
	grid_stencil(MPI_Comm comm, grid_topo&& grid);
	grid_topo & grid() { return grid_; }
	const grid_topo & grid() const { return grid_; }
	MPI_Comm comm;

private:
	grid_topo grid_;
};

}}}

#endif

