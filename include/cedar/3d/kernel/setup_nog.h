#ifndef CEDAR_3D_KERNEL_SETUP_NOG_H
#define CEDAR_3D_KERNEL_SETUP_NOG_H

#include <cedar/kernel_params.h>
#include "cedar/mpi/grid_topo.h"
#include "cedar/3d/types.h"

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	void fortran_setup_nog(const kernel_params & params, grid_topo & topo, len_t min_coarse,
		int *nog);
}

}}}
#endif
