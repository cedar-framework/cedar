#ifndef BOXMG_2D_CORE_STENCIL_OP_H
#define BOXMG_2D_CORE_STENCIL_OP_H

#include <boxmg/2d/discrete_op.h>
#include <boxmg/2d/grid_stencil.h>
#include <boxmg/types.h>
#include <boxmg/util/kernel_registry.h>

namespace boxmg { namespace bmg2d {

class stencil_op : public discrete_op
{

public:
	stencil_op() {};
stencil_op(len_t nx, len_t ny, bool intergrid) : gs(nx,ny,1,intergrid), kernels(nullptr){}
stencil_op(len_t nx, len_t ny) : gs(nx,ny), kernels(nullptr) {}
	real_t * data() { return gs.data(); }
	virtual void apply(const grid_func &x, grid_func &y) const;
	virtual void residual(const grid_func &x, const grid_func &b, grid_func &r) const;
	virtual grid_func residual(const grid_func &x, const grid_func &b) const;
	virtual void set_registry(std::shared_ptr<KernelRegistry> kreg);
	virtual std::shared_ptr<KernelRegistry> get_registry() const;
	grid_stencil & stencil() { return gs; }
	const grid_stencil & stencil() const { return gs; }
	friend std::ostream & operator<< (std::ostream &os, const stencil_op & op);

protected:
	grid_stencil gs;
	std::shared_ptr<KernelRegistry> kernels;

};

}}

#endif
