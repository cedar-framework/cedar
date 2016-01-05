#ifndef BOXMG_3D_KERNEL_RESIDUAL_H
#define BOXMG_3D_KERNEL_RESIDUAL_H

#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/stencil_op.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void residual(const stencil_op & A, const grid_func & x, const grid_func &b, grid_func & y);
}

}}}
#endif
