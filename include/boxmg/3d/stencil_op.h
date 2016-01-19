#ifndef BOXMG_3D_STENCIL_OP_H
#define BOXMG_3D_STENCIL_OP_H

#include<iostream>

#include <boxmg/types.h>
#include <boxmg/3d/stencil_op_base.h>
#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/kernel/registry.h>

namespace boxmg { namespace bmg3 {


class stencil_op : public stencil_op_base<grid_func, kernel::registry>
{
public:
	stencil_op() {};
    stencil_op(len_t nx, len_t ny, len_t nz, bool intergrid=false) : stencil_op_base(nx,ny,nz, intergrid) {}
	friend std::ostream & operator<< (std::ostream &os, const stencil_op & op);
	virtual void apply(const grid_func & x, grid_func &y) const
	{
		stencil_op_base::apply<stencil_op>(x, y);
	}
	virtual void residual(const grid_func &x, const grid_func & b, grid_func &r) const
	{
		stencil_op_base::residual<stencil_op>(x, b, r);
	}
	virtual grid_func residual(const grid_func &x, const grid_func &b) const
	{
		return stencil_op_base::residual<stencil_op>(x, b);
	}
};

}}
#endif