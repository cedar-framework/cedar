#ifndef BOXMG_2D_KERNEL_RESTRICT_H
#define BOXMG_2D_KERNEL_RESTRICT_H

#include "boxmg/2d/stencil_op.h"
#include "boxmg/2d/inter/prolong_op.h"
#include "boxmg/2d/inter/restrict_op.h"
#include "boxmg/2d/grid_func.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void fortran_restrict(const inter::restrict_op & so,
	                      const grid_func & fine,
	                      grid_func & coarse);
	void mpi_fortran_restrict(const inter::restrict_op & so,
	                          const grid_func & fine,
	                          grid_func & coarse);
}

}}}

#endif
