#ifndef BOXMG_2D_KERNEL_MATVEC_H
#define BOXMG_2D_KERNEL_MATVEC_H

#include "boxmg/2d/mpi/stencil_op.h"
#include "boxmg/2d/mpi/grid_func.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg2d::mpi;
	void matvec(const mpi::stencil_op & so,
	            const mpi::grid_func & x,
	            mpi::grid_func & y);
}

}}}
#endif
