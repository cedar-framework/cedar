#ifndef BOXMG_2D_CORE_DISCRETE_OP_H
#define BOXMG_2D_CORE_DISCRETE_OP_H

#include "grid_func.h"

namespace boxmg { namespace bmg2d { namespace core {

class DiscreteOp
{

public:
	virtual void apply(const GridFunc & x, GridFunc & y) const = 0;
	virtual void residual(const GridFunc &x, const GridFunc &b, GridFunc &r) const = 0;
	//virtual void axpy(GridFunc &
};


}}}

#endif
