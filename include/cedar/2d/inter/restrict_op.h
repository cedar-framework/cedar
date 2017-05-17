#ifndef CEDAR_2D_INTER_RESTRICT_OP_H
#define CEDAR_2D_INTER_RESTRICT_OP_H

#include <iostream>

#include "cedar/2d/stencil_op.h"
#include "cedar/2d/inter/prolong_op.h"
#include "cedar/2d/grid_func.h"


namespace cedar { namespace cdr2 { namespace inter {

class restrict_op
{

public:
restrict_op() : P(nullptr) {}
restrict_op(prolong_op * P) : P(P) {}
	void associate(prolong_op *P) { this->P = P; }

	prolong_op & getP() { return *P; }
	const prolong_op & getP() const { return *P; }

	friend std::ostream & operator<< (std::ostream &os, const restrict_op & R);

private:
	prolong_op * P;
};
}}}

#endif
