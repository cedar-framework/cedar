#ifndef BOXMG_2D_CORE_DISCRETE_OP_H
#define BOXMG_2D_CORE_DISCRETE_OP_H

#include "grid_func.h"

namespace boxmg { namespace bmg2d {

class DiscreteOp
{

public:
	virtual void apply(const grid_func & x, grid_func & y) const = 0;
	virtual void residual(const grid_func &x, const grid_func &b, grid_func &r) const = 0;
	//virtual void axpy(grid_func &
};


}}

#endif
