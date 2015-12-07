#ifndef BOXMG_2D_KERNEL_MATVEC_H
#define BOXMG_2D_KERNEL_MATVEC_H

#include "boxmg/2d/stencil_op.h"
#include "boxmg/2d/grid_func.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void matvec(const stencil_op & so,
	            const grid_func & x,
	            grid_func & y);
}

}}}
#endif
