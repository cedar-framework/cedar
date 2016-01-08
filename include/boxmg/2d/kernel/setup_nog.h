#ifndef BOXMG_2D_KERNEL_SETUP_NOG_H
#define BOXMG_2D_KERNEL_SETUP_NOG_H

#include "boxmg/mpi/grid_topo.h"
#include "boxmg/2d/types.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void fortran_setup_nog(grid_topo & topo, len_t min_coarse,
		int *nog);
}

}}}
#endif
