#ifndef BOXMG_2D_CORE_MPI_GRID_STENCIL_H
#define BOXMG_2D_CORE_MPI_GRID_STENCIL_H


#include <mpi.h>
#include "core/grid_stencil.h"
#include "core/mpi/grid_topo.h"

namespace boxmg { namespace bmg2d { namespace mpi {

class grid_stencil : public bmg2d::grid_stencil
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
