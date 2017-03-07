#ifndef BOXMG_3D_KERNEL_MATVEC_H
#define BOXMG_3D_KERNEL_MATVEC_H

#include <boxmg/3d/mpi/stencil_op.h>
#include <boxmg/3d/mpi/grid_func.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg3::mpi;
	void matvec(const mpi::stencil_op & so,
	            const mpi::grid_func & x,
	            mpi::grid_func & y);
}

}}}
#endif
