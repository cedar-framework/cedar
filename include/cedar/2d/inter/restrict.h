#ifndef CEDAR_2D_KERNEL_RESTRICT_H
#define CEDAR_2D_KERNEL_RESTRICT_H

#include <cedar/kernel_params.h>
#include "cedar/2d/stencil_op.h"
#include "cedar/2d/inter/prolong_op.h"
#include "cedar/2d/inter/restrict_op.h"
#include "cedar/2d/inter/mpi/restrict_op.h"
#include "cedar/2d/grid_func.h"
#include "cedar/2d/mpi/grid_func.h"

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	void fortran_restrict(const kernel_params & params,
	                      const inter::restrict_op & so,
	                      const grid_func & fine,
	                      grid_func & coarse);
	namespace mpi = cedar::cdr2::mpi;
	void mpi_fortran_restrict(const kernel_params & params,
	                          const inter::mpi::restrict_op & so,
	                          const mpi::grid_func & fine,
	                          mpi::grid_func & coarse);
}

}}}

#endif
