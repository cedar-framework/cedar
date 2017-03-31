#ifndef CEDAR_2D_INTER_RESTRICT_OP_H
#define CEDAR_2D_INTER_RESTRICT_OP_H

#include <iostream>

#include "cedar/2d/stencil_op.h"
#include "cedar/2d/types.h"
#include "cedar/2d/inter/prolong_op.h"
#include "cedar/2d/grid_func.h"


namespace cedar { namespace cdr2 { namespace inter {

class restrict_op : public stencil_op
{

public:
restrict_op() {}
restrict_op(prolong_op * P) : P(P) {}
	void associate(prolong_op *P) { this->P = P; }

	prolong_op & getP() { return *P; }
	const prolong_op & getP() const { return *P; }

	virtual void apply(const grid_func &x, grid_func &y) const;
	virtual void residual(const grid_func &x, const grid_func &b, grid_func &r) const{}
	friend std::ostream & operator<< (std::ostream &os, const restrict_op & R);
	friend grid_func operator* (const restrict_op & R, const grid_func &x);

private:
	prolong_op * P;
};
}}}

#endif
