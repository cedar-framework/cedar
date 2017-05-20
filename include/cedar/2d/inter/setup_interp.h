#ifndef CEDAR_2D_KERNEL_SETUP_INTERP_H
#define CEDAR_2D_KERNEL_SETUP_INTERP_H

#include <cedar/kernel_params.h>
#include "cedar/2d/grid_func.h"
#include "cedar/2d/stencil_op.h"
#include "cedar/2d/inter/prolong_op.h"
/* #include "cedar/2d/mpi/stencil_op.h" */
/* #include "cedar/2d/inter/mpi/prolong_op.h" */

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	template<class sten>
	void setup_interp(const kernel_params & params,
	                  const stencil_op<sten> & fop,
	                  const stencil_op<nine_pt> & cop,
	                  inter::prolong_op & P);

	/* namespace mpi = cedar::cdr2::mpi; */
	/* void mpi_setup_interp(const kernel_params & params, */
	/*                       int kf, int kc, int nog, */
	/*                       const mpi::stencil_op & fop, */
	/*                       const mpi::stencil_op & cop, */
	/*                       inter::mpi::prolong_op & P); */
}

}}}

#endif
