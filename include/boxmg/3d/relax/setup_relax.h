#ifndef BOXMG_3D_SETUP_RELAX_H
#define BOXMG_3D_SETUP_RELAX_H

#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/relax_stencil.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void setup_rbgs_point(const stencil_op & so,
	                      relax_stencil & sor);
}

}}}

#endif
