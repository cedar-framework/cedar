#ifndef CEDAR_3D_INTER_RESTRICT_OP_H
#define CEDAR_3D_INTER_RESTRICT_OP_H

#include <iostream>

#include <cedar/3d/inter/prolong_op.h>


namespace cedar { namespace cdr3 { namespace inter {

class restrict_op
{
public:
restrict_op() : P(nullptr) {}
    restrict_op(prolong_op *P) : P(P) {}
	void associate(prolong_op *P) { this->P = P; }
	prolong_op & getP() { return *P; }
	const prolong_op & getP() const { return *P; }

	friend std::ostream & operator<< (std::ostream &os, const restrict_op & R);
private:
	prolong_op *P;

};

}}}

#endif

