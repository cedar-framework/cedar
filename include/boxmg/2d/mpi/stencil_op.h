#ifndef BOXMG_2D_CORE_MPI_STENCIL_OP_H
#define BOXMG_2D_CORE_MPI_STENCIL_OP_H

#include <mpi.h>

#include "boxmg/2d/mpi/grid_topo.h"
#include "boxmg/2d/stencil_op.h"
#include "boxmg/2d/mpi/grid_func.h"

namespace boxmg { namespace bmg2d { namespace mpi {

class stencil_op : public bmg2d::stencil_op
{
public:
	using bmg2d::stencil_op::residual;
	stencil_op() {}
	stencil_op(topo_ptr grid);
	MPI_Comm comm;
	topo_ptr grid_ptr() { return grid_; }
	grid_topo & grid() { return *grid_; }
	const grid_topo & grid() const { return *grid_; }
	void *halo_ctx;
	virtual grid_func residual(const grid_func &x, const grid_func &b) const;
	virtual void residual(const grid_func &x, const grid_func &b, grid_func &r) const;
	friend std::ostream & operator<< (std::ostream &os, const stencil_op & op);

protected:
	topo_ptr grid_;
};

}}}

#endif
