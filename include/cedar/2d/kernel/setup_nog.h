#ifndef CEDAR_2D_KERNEL_SETUP_NOG_H
#define CEDAR_2D_KERNEL_SETUP_NOG_H

#include <cedar/kernel_params.h>
#include <cedar/mpi/grid_topo.h>

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	void fortran_setup_nog(const kernel_params & params, grid_topo & topo, len_t min_coarse,
		int *nog);
}

}}}
#endif
