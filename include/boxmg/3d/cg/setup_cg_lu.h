#ifndef BOXMG_3D_SETUP_CG_LU_H
#define BOXMG_3D_SETUP_CG_LU_H

#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/grid_func.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void setup_cg_lu(const stencil_op &so,
	                 grid_func &ABD);

}

}}}

#endif
