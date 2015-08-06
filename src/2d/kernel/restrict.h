#ifndef BOXMG_2D_KERNEL_RESTRICT_H
#define BOXMG_2D_KERNEL_RESTRICT_H

#include "core/stencil_op.h"
#include "inter/prolong_op.h"
#include "inter/restrict_op.h"
#include "core/grid_func.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void fortran_restrict(const inter::RestrictOp & so,
	                      const core::GridFunc & fine,
	                      core::GridFunc & coarse);
	void mpi_fortran_restrict(const inter::RestrictOp & so,
	                          const core::GridFunc & fine,
	                          core::GridFunc & coarse);
}

}}}

#endif
