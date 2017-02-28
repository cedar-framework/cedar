#ifndef BOXMG_2D_KERNEL_SETUP_CG_BOXMG_H
#define BOXMG_2D_KERNEL_SETUP_CG_BOXMG_H

#include <memory>

#include "boxmg/2d/mpi/stencil_op.h"
#include "boxmg/2d/solver.h"


namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg2d::mpi;
	void setup_cg_boxmg(const mpi::stencil_op & so,
	                    std::shared_ptr<config::reader> conf,
	                    std::shared_ptr<solver> *slv);
}

}}}

#endif
