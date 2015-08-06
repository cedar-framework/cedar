#include "grid_stencil.h"


using namespace boxmg::bmg2d::core::mpi;

GridStencil::GridStencil(MPI_Comm comm, GridTopo&& grd) :
	core::GridStencil(grd.ie[0] - grd.is[0] +1, grd.ie[1] - grd.is[1] +1), comm(comm), grid_(std::move(grd)) {}
