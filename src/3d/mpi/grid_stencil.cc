#include <boxmg/3d/mpi/grid_stencil.h>

using namespace boxmg::bmg3::mpi;

grid_stencil::grid_stencil(MPI_Comm comm, grid_topo&& grd) :
	boxmg::bmg3::grid_stencil(grd.ie[0] - grd.is[0] + 1,
	                          grd.ie[1] - grd.is[1] + 1,
	                          grd.ie[2] - grd.is[2] + 1),
	comm(comm),
	grid_(std::move(grd))
{}
