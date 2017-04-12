#ifndef CEDAR_2D_KERNEL_SETUP_RELAX_H
#define CEDAR_2D_KERNEL_SETUP_RELAX_H

#include <cedar/kernel_params.h>
#include "cedar/2d/stencil_op.h"
#include "cedar/2d/relax_stencil.h"
#include "cedar/2d/mpi/stencil_op.h"

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr2::mpi;
	void setup_rbgs_point(const kernel_params & params,
	                      const stencil_op & so,
	                      relax_stencil & sor);
	void setup_rbgs_x(const kernel_params & params,
	                  const stencil_op & so,
	                  relax_stencil & sor);
	void setup_rbgs_y(const kernel_params & params,
	                  const stencil_op & so,
	                  relax_stencil & sor);
	void mpi_setup_rbgs_point(const kernel_params & params,
	                          const mpi::stencil_op & so,
	                          relax_stencil & sor);
	void mpi_setup_rbgs_x(const kernel_params & params,
	                      const mpi::stencil_op & so,
	                      relax_stencil & sor);
	void mpi_setup_rbgs_y(const kernel_params & params,
	                      const mpi::stencil_op & so,
	                      relax_stencil & sor);
}

}}}

#endif
