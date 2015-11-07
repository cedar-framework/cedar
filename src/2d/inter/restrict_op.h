#ifndef BOXMG_2D_INTER_RESTRICT_OP_H
#define BOXMG_2D_INTER_RESTRICT_OP_H

#include <iostream>

#include "boxmg-common.h"

#include "core/discrete_op.h"
#include "core/types.h"
#include "prolong_op.h"
#include "core/grid_func.h"


namespace boxmg { namespace bmg2d { namespace inter {

class RestrictOp : public DiscreteOp
{

public:
RestrictOp(): kernels(nullptr) {}
RestrictOp(ProlongOp * P) : kernels(nullptr), P(P) {}
	void associate(ProlongOp *P) { this->P = P; }

	ProlongOp & getP() { return *P; }
	const ProlongOp & getP() const { return *P; }

	virtual void apply(const GridFunc &x, GridFunc &y) const;
	virtual void residual(const GridFunc &x, const GridFunc &b, GridFunc &r) const{}
	virtual void set_registry(std::shared_ptr<KernelRegistry> kreg);
	friend std::ostream & operator<< (std::ostream &os, const RestrictOp & R);
	friend GridFunc operator* (const RestrictOp & R, const GridFunc &x);

protected:
	std::shared_ptr<KernelRegistry> kernels;

private:
	ProlongOp * P;
};
}}}

#endif
