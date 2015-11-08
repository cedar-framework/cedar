#ifndef BOXMG_2D_INTER_MPI_PROLONG_OP_H
#define BOXMG_2D_INTER_MPI_PROLONG_OP_H

#include "inter/prolong_op.h"
#include "core/mpi/grid_topo.h"


namespace boxmg { namespace bmg2d { namespace inter { namespace mpi {

class prolong_op : public inter::prolong_op
{
public:
	prolong_op() {};
	prolong_op(bmg2d::mpi::topo_ptr grid);
	void *halo_ctx;
	bmg2d::mpi::grid_topo & grid() { return *grid_; }
	const bmg2d::mpi::grid_topo & grid() const { return *grid_; }
	bmg2d::mpi::topo_ptr grid_ptr() const { return grid_; }
	friend std::ostream & operator<< (std::ostream &os, const prolong_op & P);

private:
	bmg2d::mpi::topo_ptr grid_;
};

}}}}


#endif
