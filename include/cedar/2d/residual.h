#ifndef CEDAR_2D_KERNEL_RESIDUAL_H
#define CEDAR_2D_KERNEL_RESIDUAL_H

#include <cedar/kernel_params.h>
#include "cedar/2d/grid_func.h"
#include "cedar/2d/stencil_op.h"
#include "cedar/2d/mpi/stencil_op.h"
#include "cedar/2d/mpi/grid_func.h"

namespace cedar { namespace cdr2 { namespace kernel {
namespace impls
{
	namespace mpi = cedar::cdr2::mpi;
	void residual(const kernel_params & params, const stencil_op & A, const grid_func & x,
				  const grid_func & b, grid_func &y);
	void residual_fortran(const kernel_params & params, const stencil_op & A, const grid_func & x,
						  const grid_func & b, grid_func &y);
	void mpi_residual_fortran(const kernel_params & params, const mpi::stencil_op & A, const mpi::grid_func & x,
	                          const mpi::grid_func & b, mpi::grid_func &y);
}

}}}

#endif
