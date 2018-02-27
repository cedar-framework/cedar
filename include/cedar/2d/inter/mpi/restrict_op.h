#ifndef CEDAR_2D_INTER_MPI_RESTRICT_OP_H
#define CEDAR_2D_INTER_MPI_RESTRICT_OP_H

#include "cedar/2d/mpi/stencil_op.h"
#include "cedar/2d/inter/mpi/prolong_op.h"
#include "cedar/2d/mpi/grid_func.h"

namespace cedar { namespace cdr2 { namespace inter { namespace mpi {
namespace mpi = cedar::cdr2::mpi;
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

}}}}

#endif
