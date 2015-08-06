#ifndef BOXMG_2D_KERNEL_RESIDUAL_H
#define BOXMG_2D_KERNEL_RESIDUAL_H


#include "core/grid_func.h"
#include "core/stencil_op.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void residual(const core::StencilOp & A, const core::GridFunc & x,
				  const core::GridFunc & b, core::GridFunc &y);
	void residual_fortran(const core::StencilOp & A, const core::GridFunc & x,
						  const core::GridFunc & b, core::GridFunc &y);
	void mpi_residual_fortran(const core::StencilOp & A, const core::GridFunc & x,
	                          const core::GridFunc & b, core::GridFunc &y);
}

}}}

#endif
