#ifndef CEDAR_3D_SETUP_RELAX_H
#define CEDAR_3D_SETUP_RELAX_H

#include <cedar/kernel_params.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/relax_stencil.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/relax/setup_planes.h>


namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;
	void setup_rbgs_point(const kernel_params & params,
	                      const stencil_op & so,
	                      relax_stencil & sor);
	void mpi_setup_rbgs_point(const kernel_params & params,
	                          const mpi::stencil_op & so,
	                          relax_stencil & sor);
}

}}}

#endif
