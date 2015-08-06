#include "mpi_grid.h"

namespace boxmg { namespace bmg2d { namespace util { namespace mpi {

bool has_boundary(core::mpi::GridTopo & grid,core::Dir dir)
{
	if (dir == core::Dir::N)
		return grid.is(1) - 1 + grid.nlocal(1) == grid.nglobal(1);
	if (dir == core::Dir::S)
		return grid.is(1) == 1;
	if (dir == core::Dir::E)
		return grid.is(0) - 1 + grid.nlocal(0) == grid.nglobal(0);
	if (dir == core::Dir::W)
		return grid.is(0) == 1;
	else {
		log::error << "Invalid direction (boundary check)" << std::endl;
		return false;
	}
}

}}}}
