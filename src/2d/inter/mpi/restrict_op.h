#ifndef BOXMG_2D_INTER_MPI_RESTRICT_OP_H
#define BOXMG_2D_INTER_MPI_RESTRICT_OP_H

#include "inter/restrict_op.h"
#include "core/mpi/grid_func.h"

namespace boxmg { namespace bmg2d { namespace inter { namespace mpi {

class RestrictOp : public inter::RestrictOp
{
public:
	friend bmg2d::mpi::GridFunc operator*(const RestrictOp &R, const GridFunc &x);
};

}}}}

#endif
