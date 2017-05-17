#ifndef CEDAR_LEVEL_H
#define CEDAR_LEVEL_H

#include <memory>


namespace cedar {
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

}

#endif
