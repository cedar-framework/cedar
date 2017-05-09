#ifndef CEDAR_2D_KERNEL_SETUP_INTERP_H
#define CEDAR_2D_KERNEL_SETUP_INTERP_H

#include <cedar/kernel_params.h>
#include "cedar/2d/grid_func.h"
#include "cedar/2d/stencil_op.h"
#include "cedar/2d/inter/prolong_op.h"
#include "cedar/2d/mpi/stencil_op.h"
#include "cedar/2d/inter/mpi/prolong_op.h"

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr2::mpi;
	void setup_interp(const kernel_params & params,
	                  int kf, int kc, int nog,
	                  const stencil_op & fop,
	                  const stencil_op & cop,
	                  inter::prolong_op & P);
	void mpi_setup_interp(const kernel_params & params,
	                      int kf, int kc, int nog,
	                      const mpi::stencil_op & fop,
	                      const mpi::stencil_op & cop,
	                      inter::mpi::prolong_op & P);
}

}}}

#endif
