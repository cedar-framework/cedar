#ifndef CEDAR_2D_KERNEL_SETUP_CG_LU_H
#define CEDAR_2D_KERNEL_SETUP_CG_LU_H

#include <cedar/kernel_params.h>
#include "cedar/2d/stencil_op.h"
#include "cedar/2d/grid_func.h"
//#include "cedar/2d/mpi/stencil_op.h"

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	template<class sten>
	void setup_cg_lu(const kernel_params & params, const stencil_op<sten> & so,
	                 grid_func & ABD);
	/* namespace mpi = cedar::cdr2::mpi; */
	/* void mpi_setup_cg_lu(const kernel_params & params, const mpi::stencil_op & so, */
	/*                      grid_func & ABD); */
}

}}}

#endif
