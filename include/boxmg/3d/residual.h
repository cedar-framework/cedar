#ifndef BOXMG_3D_KERNEL_RESIDUAL_H
#define BOXMG_3D_KERNEL_RESIDUAL_H

#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/mpi/stencil_op.h>
#include <boxmg/3d/mpi/grid_func.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg3::mpi;
	void residual(const stencil_op & A, const grid_func & x, const grid_func &b, grid_func & y);

	void mpi_residual_fortran(const mpi::stencil_op & A, const mpi::grid_func & x,
	                          const mpi::grid_func & b, mpi::grid_func &y);
}

}}}
#endif
