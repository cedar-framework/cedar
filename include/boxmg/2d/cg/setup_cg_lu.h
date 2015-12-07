#ifndef BOXMG_2D_KERNEL_SETUP_CG_LU_H
#define BOXMG_2D_KERNEL_SETUP_CG_LU_H

#include "boxmg/2d/stencil_op.h"
#include "boxmg/2d/grid_func.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void setup_cg_lu(const stencil_op & so,
	                 grid_func & ABD);

	void mpi_setup_cg_lu(const stencil_op & so,
	                     grid_func & ABD);

}

}}}

#endif
