#ifndef CEDAR_3D_KERNEL_RESIDUAL_H
#define CEDAR_3D_KERNEL_RESIDUAL_H

#include <cedar/3d/grid_func.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/mpi/grid_func.h>

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;
	void residual(const stencil_op & A, const grid_func & x, const grid_func &b, grid_func & y);

	void mpi_residual_fortran(const mpi::stencil_op & A, const mpi::grid_func & x,
	                          const mpi::grid_func & b, mpi::grid_func &y);
}

}}}
#endif
