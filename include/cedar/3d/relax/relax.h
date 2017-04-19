#ifndef CEDAR_3D_RELAX_H
#define CEDAR_3D_RELAX_H

#include <cedar/kernel_params.h>
#include <cedar/cycle/types.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/grid_func.h>
#include <cedar/3d/relax_stencil.h>
#include <cedar/3d/mpi/grid_func.h>
#include <cedar/3d/mpi/stencil_op.h>

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;
	void relax_rbgs_point(const kernel_params & params,
	                      const stencil_op & so,
	                      grid_func &x,
	                      const grid_func &b,
	                      const relax_stencil & sor,
	                      cycle::Dir cycle_dir);

	void mpi_relax_rbgs_point(const kernel_params & params,
	                          const mpi::stencil_op & so,
	                          mpi::grid_func & x,
	                          const mpi::grid_func & b,
	                          const relax_stencil & sor,
	                          cycle::Dir cycle_dir);
}

}}}

#endif
