#ifndef BOXMG_2D_KERNEL_SETUP_CG_LU_H
#define BOXMG_2D_KERNEL_SETUP_CG_LU_H

#include "boxmg/2d/stencil_op.h"
#include "boxmg/2d/grid_func.h"
#include "boxmg/2d/mpi/stencil_op.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg2d::mpi;
	void setup_cg_lu(const stencil_op & so,
	                 grid_func & ABD);

	void mpi_setup_cg_lu(const mpi::stencil_op & so,
	                     grid_func & ABD);

}

}}}

#endif
