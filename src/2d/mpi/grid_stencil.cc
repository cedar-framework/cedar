#include <boxmg/2d/mpi/grid_stencil.h>


using namespace boxmg::bmg2d::mpi;

grid_stencil::grid_stencil(MPI_Comm comm, grid_topo&& grd) :
	boxmg::bmg2d::grid_stencil(grd.ie[0] - grd.is[0] +1, grd.ie[1] - grd.is[1] +1), comm(comm), grid_(std::move(grd)) {}
