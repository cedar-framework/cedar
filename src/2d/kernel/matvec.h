#ifndef BOXMG_2D_KERNEL_MATVEC_H
#define BOXMG_2D_KERNEL_MATVEC_H

#include "core/stencil_op.h"
#include "core/grid_func.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void matvec(const core::StencilOp & so,
	            const core::GridFunc & x,
	            core::GridFunc y);
}

}}}
#endif
