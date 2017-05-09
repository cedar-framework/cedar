#ifndef CEDAR_3D_RESTRICT_H
#define CEDAR_3D_RESTRICT_H

#include <cedar/kernel_params.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/inter/prolong_op.h>
#include <cedar/3d/inter/restrict_op.h>
#include <cedar/3d/grid_func.h>
#include <cedar/3d/inter/mpi/restrict_op.h>
#include <cedar/3d/mpi/grid_func.h>

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;

	void fortran_restrict(const kernel_params & params,
	                      const inter::restrict_op & so,
	                      const grid_func & fine,
	                      grid_func & coarse);
	void mpi_fortran_restrict(const kernel_params & params,
	                          const inter::mpi::restrict_op & so,
	                          const mpi::grid_func & fine,
	                          mpi::grid_func & coarse);
}

}}}

#endif
