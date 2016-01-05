#ifndef BOXMG_3D_INTER_SETUP_INTERP_H
#define BOXMG_3D_INTER_SETUP_INTERP_H

#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/inter/prolong_op.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void setup_interp(int kf, int kc, int nog,
	                  const stencil_op & fop,
	                  const stencil_op & cop,
	                  inter::prolong_op & P);
}

}}}

#endif
