#ifndef BOXMG_2D_KERNEL_SETUP_NOG_H
#define BOXMG_2D_KERNEL_SETUP_NOG_H

#include "core/mpi/grid_topo.h"
#include "core/types.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void fortran_setup_nog(bmg2d::mpi::GridTopo & topo, len_t min_coarse,
		int *nog);
}

}}}
#endif
