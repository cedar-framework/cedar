#ifndef BOXMG_2D_KERNEL_SETUP_INTERP_H
#define BOXMG_2D_KERNEL_SETUP_INTERP_H

#include "core/grid_func.h"
#include "core/stencil_op.h"
#include "inter/prolong_op.h"
#include "core/mpi/stencil_op.h"
#include "inter/mpi/prolong_op.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void setup_interp(int kf, int kc, int nog,
	                  const StencilOp & fop,
	                  const StencilOp & cop,
	                  inter::prolong_op & P);
	void mpi_setup_interp(int kf, int kc, int nog,
	                      const StencilOp & fop,
	                      const StencilOp & cop,
	                      inter::prolong_op & P);
}

}}}

#endif
