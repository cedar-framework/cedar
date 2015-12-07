#ifndef BOXMG_2D_INTER_RESTRICT_OP_H
#define BOXMG_2D_INTER_RESTRICT_OP_H

#include <iostream>

#include "boxmg/util/kernel_registry.h"

#include "boxmg/2d/discrete_op.h"
#include "boxmg/2d/types.h"
#include "boxmg/2d/inter/prolong_op.h"
#include "boxmg/2d/grid_func.h"


namespace boxmg { namespace bmg2d { namespace inter {

class restrict_op : public discrete_op
{

public:
restrict_op(): kernels(nullptr) {}
restrict_op(prolong_op * P) : kernels(nullptr), P(P) {}
	void associate(prolong_op *P) { this->P = P; }

	prolong_op & getP() { return *P; }
	const prolong_op & getP() const { return *P; }

	virtual void apply(const grid_func &x, grid_func &y) const;
	virtual void residual(const grid_func &x, const grid_func &b, grid_func &r) const{}
	virtual void set_registry(std::shared_ptr<KernelRegistry> kreg);
	friend std::ostream & operator<< (std::ostream &os, const restrict_op & R);
	friend grid_func operator* (const restrict_op & R, const grid_func &x);

protected:
	std::shared_ptr<KernelRegistry> kernels;

private:
	prolong_op * P;
};
}}}

#endif
