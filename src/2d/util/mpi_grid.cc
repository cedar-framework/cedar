#include "mpi_grid.h"

namespace boxmg { namespace bmg2d { namespace util { namespace mpi {

bool has_boundary(bmg2d::mpi::grid_topo & grid,bmg2d::Dir dir)
{
	if (dir == bmg2d::Dir::N)
		return grid.is(1) - 1 + grid.nlocal(1) == grid.nglobal(1);
	if (dir == bmg2d::Dir::S)
		return grid.is(1) == 1;
	if (dir == bmg2d::Dir::E)
		return grid.is(0) - 1 + grid.nlocal(0) == grid.nglobal(0);
	if (dir == bmg2d::Dir::W)
		return grid.is(0) == 1;
	else {
		log::error << "Invalid direction (boundary check)" << std::endl;
		return false;
	}
}

}}}}
