#ifndef BOXMG_2D_CORE_MPI_STENCIL_OP_H
#define BOXMG_2D_CORE_MPI_STENCIL_OP_H

#include <mpi.h>

#include <boxmg/2d/stencil_op_base.h>
#include <boxmg/mpi/grid_topo.h>
#include <boxmg/2d/mpi/grid_func.h>
#include <boxmg/2d/kernel/mpi/registry.h>

namespace boxmg { namespace bmg2d { namespace mpi {

class stencil_op : public stencil_op_base<grid_func, kernel::mpi::registry>
{
public:
	stencil_op() {}
	stencil_op(topo_ptr grid);
    stencil_op(len_t nx, len_t ny, bool intergrid=false): stencil_op_base(nx,ny,intergrid) {}
	MPI_Comm comm;
	topo_ptr grid_ptr() { return grid_; }
	grid_topo & grid() { return *grid_; }
	const grid_topo & grid() const { return *grid_; }
	void *halo_ctx;
	virtual grid_func residual(const grid_func &x, const grid_func &b) const;
	friend std::ostream & operator<< (std::ostream &os, const stencil_op & op);
	virtual void apply(const grid_func &x, grid_func &y) const
	{
		stencil_op_base::apply<stencil_op>(x, y);
	}
	virtual void residual(const grid_func &x, const grid_func & b, grid_func &r) const
	{
		stencil_op_base::residual<stencil_op>(x, b, r);
	}

protected:
	topo_ptr grid_;
};

}}}

#endif
