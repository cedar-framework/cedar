#ifndef BOXMG_2D_KERNEL_SETUP_INTERP_H
#define BOXMG_2D_KERNEL_SETUP_INTERP_H

#include "boxmg/2d/grid_func.h"
#include "boxmg/2d/stencil_op.h"
#include "boxmg/2d/inter/prolong_op.h"
#include "boxmg/2d/mpi/stencil_op.h"
#include "boxmg/2d/inter/mpi/prolong_op.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg2d::mpi;
	void setup_interp(int kf, int kc, int nog,
	                  const stencil_op & fop,
	                  const stencil_op & cop,
	                  inter::prolong_op & P);
	void mpi_setup_interp(int kf, int kc, int nog,
	                      const mpi::stencil_op & fop,
	                      const mpi::stencil_op & cop,
	                      inter::mpi::prolong_op & P);
}

}}}

#endif
