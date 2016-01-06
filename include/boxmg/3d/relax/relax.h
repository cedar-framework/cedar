#ifndef BOXMG_3D_RELAX_H
#define BOXMG_3D_RELAX_H

#include <boxmg/cycle/types.h>
#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/relax_stencil.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void relax_rbgs_point(const stencil_op & so,
	                      grid_func &x,
	                      const grid_func &b,
	                      const relax_stencil & sor,
	                      cycle::Dir cycle_dir);
}

}}}

#endif
