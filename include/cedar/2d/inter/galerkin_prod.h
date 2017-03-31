#ifndef CEDAR_2D_KERNEL_GALERKIN_PROD_H
#define CEDAR_2D_KERNEL_GALERKIN_PROD_H

#include "cedar/2d/stencil_op.h"
#include "cedar/2d/inter/prolong_op.h"
#include "cedar/2d/grid_func.h"
#include "cedar/2d/mpi/stencil_op.h"
#include "cedar/2d/inter/mpi/prolong_op.h"

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr2::mpi;
	void galerkin_prod(int kf, int kc, int nog,
	                   const inter::prolong_op & P,
	                   const stencil_op & fop,
	                   stencil_op & cop);


	void mpi_galerkin_prod(int kf, int kc, int nog,
	                       const inter::mpi::prolong_op & P,
	                       const mpi::stencil_op & fop,
	                       mpi::stencil_op & cop);
}

}}}

#endif