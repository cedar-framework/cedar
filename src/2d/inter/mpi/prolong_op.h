#ifndef BOXMG_2D_INTER_MPI_PROLONG_OP_H
#define BOXMG_2D_INTER_MPI_PROLONG_OP_H

#include "inter/prolong_op.h"
#include "core/mpi/grid_topo.h"


namespace boxmg { namespace bmg2d { namespace inter { namespace mpi {

class ProlongOp : public inter::ProlongOp
{
public:
	ProlongOp() {};
	ProlongOp(bmg2d::mpi::topo_ptr grid);
	void *halo_ctx;
	bmg2d::mpi::GridTopo & grid() { return *grid_; }
	const bmg2d::mpi::GridTopo & grid() const { return *grid_; }
	bmg2d::mpi::topo_ptr grid_ptr() const { return grid_; }
	friend std::ostream & operator<< (std::ostream &os, const ProlongOp & P);

private:
	bmg2d::mpi::topo_ptr grid_;
};

}}}}


#endif
