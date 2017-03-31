#ifndef CEDAR_3D_SETUP_CG_LU_H
#define CEDAR_3D_SETUP_CG_LU_H

#include <cedar/3d/stencil_op.h>
#include <cedar/3d/grid_func.h>
#include <cedar/3d/mpi/stencil_op.h>

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;
	void setup_cg_lu(const stencil_op &so,
	                 grid_func &ABD);

	void mpi_setup_cg_lu(const mpi::stencil_op & so,
	                     grid_func & ABD);
}

}}}

#endif
