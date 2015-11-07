#ifndef BOXMG_2D_SOLVER_LEVEL_H
#define BOXMG_2D_SOLVER_LEVEL_H

#include <memory>

#include "core/discrete_op.h"

namespace boxmg { namespace bmg2d {

struct Level
{
	Level() {}
    /* Level(core::DiscreteOp &A, core::DiscreteOp &P) :A(A),P(P) {} */
	/* core::DiscreteOp & A; */
	/* core::DiscreteOp & P; */
	std::function<void(const DiscreteOp & A, GridFunc &x, const GridFunc &b)> presmoother;
	std::function<void(const DiscreteOp & A, GridFunc &x, const GridFunc &b)> postsmoother;
	//DiscreteOp & R;
};

}}

#endif
