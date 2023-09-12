#include <cedar/2d/gpu/grid_stencil.h>


using namespace cedar::cdr2::gpu::mpi;

grid_stencil::grid_stencil(MPI_Comm comm, grid_topo&& grd) :
	cedar::cdr2::grid_stencil(grd.ie[0] - grd.is[0] +1, grd.ie[1] - grd.is[1] +1), comm(comm), grid_(std::move(grd)) {}
