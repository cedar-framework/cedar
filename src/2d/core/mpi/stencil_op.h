#ifndef BOXMG_2D_CORE_MPI_STENCIL_OP_H
#define BOXMG_2D_CORE_MPI_STENCIL_OP_H

#include <mpi.h>

#include "grid_topo.h"
#include "core/stencil_op.h"
#include "core/mpi/grid_func.h"

namespace boxmg { namespace bmg2d { namespace mpi {

class StencilOp : public bmg2d::StencilOp
{
public:
	using bmg2d::StencilOp::residual;
	StencilOp() {}
	StencilOp(topo_ptr grid);
	MPI_Comm comm;
	topo_ptr grid_ptr() { return grid_; }
	GridTopo & grid() { return *grid_; }
	const GridTopo & grid() const { return *grid_; }
	void *halo_ctx;
	virtual grid_func residual(const grid_func &x, const grid_func &b) const;
	virtual void residual(const grid_func &x, const grid_func &b, grid_func &r) const;
	friend std::ostream & operator<< (std::ostream &os, const StencilOp & op);

protected:
	topo_ptr grid_;
};

}}}

#endif
