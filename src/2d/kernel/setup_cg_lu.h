#ifndef BOXMG_2D_KERNEL_SETUP_CG_LU_H
#define BOXMG_2D_KERNEL_SETUP_CG_LU_H

#include "core/stencil_op.h"
#include "core/grid_func.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void setup_cg_lu(const StencilOp & so,
	                 grid_func & ABD);

	void mpi_setup_cg_lu(const StencilOp & so,
	                     grid_func & ABD);

}

}}}

#endif
