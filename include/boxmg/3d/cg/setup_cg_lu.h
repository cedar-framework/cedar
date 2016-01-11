#ifndef BOXMG_3D_SETUP_CG_LU_H
#define BOXMG_3D_SETUP_CG_LU_H

#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/mpi/stencil_op.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg3::mpi;
	void setup_cg_lu(const stencil_op &so,
	                 grid_func &ABD);

	void mpi_setup_cg_lu(const mpi::stencil_op & so,
	                     grid_func & ABD);
}

}}}

#endif
