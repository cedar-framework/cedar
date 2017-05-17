#ifndef CEDAR_2D_CORE_MPI_STENCIL_OP_H
#define CEDAR_2D_CORE_MPI_STENCIL_OP_H

#include <mpi.h>

#include <cedar/mpi/par_object.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/kernel/mpi/registry.h>

namespace cedar { namespace cdr2 { namespace mpi {

class stencil_op : public stencil_op_base<grid_func, kernel::mpi::registry>, public par_object
{
public:
	stencil_op() {}
	stencil_op(topo_ptr grid);
    stencil_op(len_t nx, len_t ny, bool intergrid=false): stencil_op_base(nx,ny,intergrid) {}
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
};

}}}

#endif
