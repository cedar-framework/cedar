#ifndef CEDAR_2D_KERNEL_MATVEC_H
#define CEDAR_2D_KERNEL_MATVEC_H

#include "cedar/2d/mpi/stencil_op.h"
#include "cedar/2d/mpi/grid_func.h"

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr2::mpi;
	void matvec(const mpi::stencil_op & so,
	            const mpi::grid_func & x,
	            mpi::grid_func & y);
}

}}}
#endif
