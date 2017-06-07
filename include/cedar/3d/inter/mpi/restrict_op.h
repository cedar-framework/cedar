#ifndef CEDAR_3D_INTER_MPI_RESTRICT_OP_H
#define CEDAR_3D_INTER_MPI_RESTRICT_OP_H

#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/inter/mpi/prolong_op.h>
#include <cedar/3d/mpi/grid_func.h>

namespace cedar { namespace cdr3 { namespace inter { namespace mpi {
namespace mpi = cedar::cdr3::mpi;
class restrict_op
{
public:
restrict_op() {}
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
