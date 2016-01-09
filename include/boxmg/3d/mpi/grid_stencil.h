#ifndef BOXMG_3D_MPI_GRID_STENCIL_H
#define BOXMG_3D_MPI_GRID_STENCIL_H

#include <mpi.h>
#include <boxmg/3d/grid_stencil.h>
#include <boxmg/mpi/grid_topo.h>

namespace boxmg { namespace bmg3 { namespace mpi {

class grid_stencil : public bmg3::grid_stencil
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

