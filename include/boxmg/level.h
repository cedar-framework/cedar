#ifndef BOXMG_LEVEL_H
#define BOXMG_LEVEL_H

#include <memory>

#include "boxmg/discrete_op.h"

namespace boxmg { namespace bmg2d {
template <class grid_func>
struct Level
{
	Level() {}
    /* Level(core::discrete_op &A, core::discrete_op &P) :A(A),P(P) {} */
	/* core::discrete_op & A; */
	/* core::discrete_op & P; */
	std::function<void(const discrete_op<grid_func> & A, grid_func &x, const grid_func &b)> presmoother;
	std::function<void(const discrete_op<grid_func> & A, grid_func &x, const grid_func &b)> postsmoother;
	//discrete_op & R;
};

}}

#endif
