#ifndef CEDAR_DISCRETE_OP_H
#define CEDAR_DISCRETE_OP_H


namespace cedar {

template <class grid_func>
class discrete_op
{

public:
	virtual void apply(const grid_func & x, grid_func & y) const = 0;
	virtual void residual(const grid_func &x, const grid_func &b, grid_func &r) const = 0;
	//virtual void axpy(grid_func &
};

}

#endif
