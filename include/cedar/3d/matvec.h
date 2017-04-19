#ifndef CEDAR_3D_KERNEL_MATVEC_H
#define CEDAR_3D_KERNEL_MATVEC_H

#include <cedar/kernel_params.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/mpi/grid_func.h>

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;
	void matvec(const kernel_params & params,
	            const mpi::stencil_op & so,
	            const mpi::grid_func & x,
	            mpi::grid_func & y);
}

}}}
#endif
