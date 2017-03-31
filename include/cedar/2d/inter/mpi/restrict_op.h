#ifndef CEDAR_2D_INTER_MPI_RESTRICT_OP_H
#define CEDAR_2D_INTER_MPI_RESTRICT_OP_H

#include "cedar/2d/mpi/stencil_op.h"
#include "cedar/2d/inter/mpi/prolong_op.h"
#include "cedar/2d/mpi/grid_func.h"

namespace cedar { namespace cdr2 { namespace inter { namespace mpi {
namespace mpi = cedar::cdr2::mpi;
class restrict_op : public mpi::stencil_op
{
public:
restrict_op() {}
restrict_op(prolong_op * P) : P(P) {}
	void associate(prolong_op *P) { this->P = P; }
	prolong_op & getP() { return *P; }
	const prolong_op & getP() const { return *P; }
	virtual void apply(const mpi::grid_func &x, mpi::grid_func &y) const;
	virtual void residual(const mpi::grid_func &x, const mpi::grid_func &b, mpi::grid_func &r) const{}
	friend mpi::grid_func operator*(const restrict_op &R, const mpi::grid_func &x);
private:
	prolong_op * P;
};

}}}}

#endif
