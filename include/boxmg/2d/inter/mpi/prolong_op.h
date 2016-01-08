#ifndef BOXMG_2D_INTER_MPI_PROLONG_OP_H
#define BOXMG_2D_INTER_MPI_PROLONG_OP_H

#include "boxmg/2d/mpi/stencil_op.h"
#include "boxmg/2d/mpi/grid_func.h"


namespace boxmg { namespace bmg2d { namespace inter { namespace mpi {


namespace mpi = boxmg::bmg2d::mpi;
using iadd_pack = std::tuple<const prolong_op&, const mpi::grid_func &, const mpi::grid_func&>;

class prolong_op : public mpi::stencil_op
{
public:
	prolong_op() {};
	prolong_op(topo_ptr grid);
	friend std::ostream & operator<< (std::ostream &os, const prolong_op & P);
	mpi::stencil_op * fine_op;
	mpi::grid_func *residual;
	friend iadd_pack operator*(const prolong_op&, const mpi::grid_func&);
};

}}}}


#endif
