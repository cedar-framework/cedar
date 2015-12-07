#ifndef BOXMG_2D_KERNEL_INTERP_H
#define BOXMG_2D_KERNEL_INTERP_H

#include "boxmg/2d/stencil_op.h"
#include "boxmg/2d/inter/prolong_op.h"
#include "boxmg/2d/grid_func.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void fortran_interp(const inter::prolong_op & P,
	                    const grid_func & coarse,
	                    const grid_func & residual,
	                    grid_func & fine);
	void mpi_fortran_interp(const inter::prolong_op & P,
	                        const grid_func & coarse,
	                        const grid_func & residual,
	                        grid_func & fine);
}

}}}

#endif
