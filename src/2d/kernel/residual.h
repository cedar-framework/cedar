#ifndef BOXMG_2D_KERNEL_RESIDUAL_H
#define BOXMG_2D_KERNEL_RESIDUAL_H


#include "core/grid_func.h"
#include "core/stencil_op.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void residual(const StencilOp & A, const grid_func & x,
				  const grid_func & b, grid_func &y);
	void residual_fortran(const StencilOp & A, const grid_func & x,
						  const grid_func & b, grid_func &y);
	void mpi_residual_fortran(const StencilOp & A, const grid_func & x,
	                          const grid_func & b, grid_func &y);
}

}}}

#endif
