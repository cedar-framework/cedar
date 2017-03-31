#include "cedar/2d/util/mpi_grid.h"

namespace cedar { namespace cdr2 { namespace util { namespace mpi {

bool has_boundary(grid_topo & grid,cdr2::dir dir)
{
	if (dir == cdr2::dir::N)
		return grid.is(1) - 1 + grid.nlocal(1) == grid.nglobal(1);
	if (dir == cdr2::dir::S)
		return grid.is(1) == 1;
	if (dir == cdr2::dir::E)
		return grid.is(0) - 1 + grid.nlocal(0) == grid.nglobal(0);
	if (dir == cdr2::dir::W)
		return grid.is(0) == 1;
	else {
		log::error << "Invalid direction (boundary check)" << std::endl;
		return false;
	}
}

}}}}
