#ifndef BOXMG_2D_KERNEL_SETUP_CG_BOXMG_H
#define BOXMG_2D_KERNEL_SETUP_CG_BOXMG_H

#include <memory>

#include "core/stencil_op.h"
#include "core/solver.h"


namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{

	void setup_cg_boxmg(const stencil_op & so,
	                    std::shared_ptr<solver> *slv);
}

}}}

#endif
