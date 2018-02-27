#ifndef CEDAR_2D_KERNEL_GALERKIN_PROD_H
#define CEDAR_2D_KERNEL_GALERKIN_PROD_H

#include <cedar/kernel_params.h>
#include "cedar/2d/stencil_op.h"
#include "cedar/2d/inter/prolong_op.h"
#include "cedar/2d/grid_func.h"
#include "cedar/2d/mpi/stencil_op.h"
#include "cedar/2d/inter/mpi/prolong_op.h"

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	template<class sten>
	void galerkin_prod(const kernel_params & params,
	                   const inter::prolong_op & P,
	                   const stencil_op<sten> & fop,
	                   stencil_op<nine_pt> & cop);
}

}}}

#endif
