#ifndef BOXMG_2D_SOLVER_LEVEL_H
#define BOXMG_2D_SOLVER_LEVEL_H

#include <memory>

#include "core/discrete_op.h"

namespace boxmg { namespace bmg2d { namespace solver {

struct Level
{
	Level() {}
    /* Level(core::DiscreteOp &A, core::DiscreteOp &P) :A(A),P(P) {} */
	/* core::DiscreteOp & A; */
	/* core::DiscreteOp & P; */
	std::function<void(const core::DiscreteOp & A, core::GridFunc &x, const core::GridFunc &b)> presmoother;
	std::function<void(const core::DiscreteOp & A, core::GridFunc &x, const core::GridFunc &b)> postsmoother;
	//core::DiscreteOp & R;
};

}}}

#endif
