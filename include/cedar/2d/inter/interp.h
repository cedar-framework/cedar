#ifndef CEDAR_2D_KERNEL_INTERP_H
#define CEDAR_2D_KERNEL_INTERP_H

#include <cedar/kernel_params.h>
#include "cedar/2d/stencil_op.h"
#include "cedar/2d/inter/prolong_op.h"
#include "cedar/2d/grid_func.h"
#include "cedar/2d/mpi/grid_func.h"
#include "cedar/2d/inter/mpi/prolong_op.h"

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr2::mpi;
	void fortran_interp(const kernel_params & params,
	                    const inter::prolong_op & P,
	                    const grid_func & coarse,
	                    const grid_func & residual,
	                    grid_func & fine);
	void mpi_fortran_interp(const kernel_params & params,
	                        const inter::mpi::prolong_op & P,
	                        const mpi::grid_func & coarse,
	                        const mpi::grid_func & residual,
	                        mpi::grid_func & fine);
}

}}}

#endif
