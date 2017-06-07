#ifndef CEDAR_3D_INTERP_H
#define CEDAR_3D_INTERP_H

#include <cedar/kernel_params.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/inter/prolong_op.h>
#include <cedar/3d/grid_func.h>

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	void fortran_interp(const kernel_params & params,
	                    const inter::prolong_op & P,
	                    const grid_func & coarse,
	                    const grid_func & residual,
	                    grid_func & fine);
}

}}}

#endif
