#ifndef CEDAR_2D_CORE_MPI_GRID_STENCIL_H
#define CEDAR_2D_CORE_MPI_GRID_STENCIL_H


#include <mpi.h>
#include "cedar/2d/grid_stencil.h"
#include "cedar/2d/mpi/grid_topo.h"

namespace cedar { namespace cdr2 { namespace mpi {

class grid_stencil : public cdr2::grid_stencil
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
