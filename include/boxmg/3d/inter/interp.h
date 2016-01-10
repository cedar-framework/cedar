#ifndef BOXMG_3D_INTERP_H
#define BOXMG_3D_INTERP_H

#include "boxmg/3d/stencil_op.h"
#include "boxmg/3d/inter/prolong_op.h"
#include "boxmg/3d/grid_func.h"

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void fortran_interp(const inter::prolong_op & P,
	                    const grid_func & coarse,
	                    const grid_func & residual,
	                    grid_func & fine);
}

}}}

#endif
