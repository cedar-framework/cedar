#ifndef BOXMG_2D_INTER_MPI_RESTRICT_OP_H
#define BOXMG_2D_INTER_MPI_RESTRICT_OP_H

#include "boxmg/2d/inter/restrict_op.h"
#include "boxmg/2d/mpi/grid_func.h"

namespace boxmg { namespace bmg2d { namespace inter { namespace mpi {

class restrict_op : public inter::restrict_op
{
public:
	friend bmg2d::mpi::grid_func operator*(const restrict_op &R, const grid_func &x);
};

}}}}

#endif
