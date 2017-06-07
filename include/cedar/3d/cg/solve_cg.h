#ifndef CEDAR_3D_SOLVE_CG_H
#define CEDAR_3D_SOLVE_CG_H

#include <cedar/3d/grid_func.h>
#include <cedar/kernel_params.h>

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	void fortran_solve_cg(const kernel_params & params,
	                      grid_func & x,
	                      const grid_func & b,
	                      const grid_func & ABD,
	                      real_t *bbd);
}

}}}

#endif
