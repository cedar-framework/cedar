#ifndef BOXMG_2D_CORE_MPI_GRID_STENCIL_H
#define BOXMG_2D_CORE_MPI_GRID_STENCIL_H


#include <mpi.h>
#include "core/grid_stencil.h"
#include "core/mpi/grid_topo.h"

namespace boxmg { namespace bmg2d { namespace mpi {

class GridStencil : public bmg2d::GridStencil
{
public:
	GridStencil();
	GridStencil(MPI_Comm comm, GridTopo&& grid);
	GridTopo & grid() { return grid_; }
	const GridTopo & grid() const { return grid_; }
	MPI_Comm comm;

private:
	GridTopo grid_;
};

}}}

#endif
