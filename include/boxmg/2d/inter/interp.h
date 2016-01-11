#ifndef BOXMG_2D_KERNEL_INTERP_H
#define BOXMG_2D_KERNEL_INTERP_H

#include "boxmg/2d/stencil_op.h"
#include "boxmg/2d/inter/prolong_op.h"
#include "boxmg/2d/grid_func.h"
#include "boxmg/2d/mpi/grid_func.h"
#include "boxmg/2d/inter/mpi/prolong_op.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg2d::mpi;
	void fortran_interp(const inter::prolong_op & P,
	                    const grid_func & coarse,
	                    const grid_func & residual,
	                    grid_func & fine);
	void mpi_fortran_interp(const inter::mpi::prolong_op & P,
	                        const mpi::grid_func & coarse,
	                        const mpi::grid_func & residual,
	                        mpi::grid_func & fine);
}

}}}

#endif
