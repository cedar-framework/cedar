#ifndef BOXMG_3D_SOLVE_CG_H
#define BOXMG_3D_SOLVE_CG_H

#include <boxmg/3d/grid_func.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void fortran_solve_cg(grid_func & x,
	                      const grid_func & b,
	                      const grid_func & ABD,
	                      real_t *bbd);
}

}}}

#endif
