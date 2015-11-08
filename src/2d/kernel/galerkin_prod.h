#ifndef BOXMG_2D_KERNEL_GALERKIN_PROD_H
#define BOXMG_2D_KERNEL_GALERKIN_PROD_H

#include "core/stencil_op.h"
#include "inter/prolong_op.h"
#include "core/grid_func.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void galerkin_prod(int kf, int kc, int nog,
	                   const inter::prolong_op & P,
	                   const StencilOp & fop,
	                   StencilOp & cop);


	void mpi_galerkin_prod(int kf, int kc, int nog,
	                       const inter::prolong_op & P,
	                       const StencilOp & fop,
	                       StencilOp & cop);
}

}}}

#endif
