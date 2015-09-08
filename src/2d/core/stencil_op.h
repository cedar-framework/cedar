#ifndef BOXMG_2D_CORE_STENCIL_OP_H
#define BOXMG_2D_CORE_STENCIL_OP_H

#include "discrete_op.h"
#include "grid_stencil.h"
#include "boxmg-common.h"

namespace boxmg { namespace bmg2d { namespace core {

class StencilOp : public DiscreteOp
{

public:
	StencilOp() {};
StencilOp(len_t nx, len_t ny, bool intergrid) : gs(nx,ny,1,intergrid), kernels(nullptr){}
StencilOp(len_t nx, len_t ny) : gs(nx,ny), kernels(nullptr) {}
	real_t * data() { return gs.data(); }
	virtual void apply(const GridFunc &x, GridFunc &y) const;
	virtual void residual(const GridFunc &x, const GridFunc &b, GridFunc &r) const;
	virtual GridFunc residual(const core::GridFunc &x, const core::GridFunc &b) const;
	virtual void set_registry(std::shared_ptr<KernelRegistry> kreg);
	virtual std::shared_ptr<KernelRegistry> get_registry() const;
	GridStencil & stencil() { return gs; }
	const GridStencil & stencil() const { return gs; }
	friend std::ostream & operator<< (std::ostream &os, const StencilOp & op);

protected:
	GridStencil gs;
	std::shared_ptr<KernelRegistry> kernels;

};

}}}

#endif
