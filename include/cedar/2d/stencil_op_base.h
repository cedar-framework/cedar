#ifndef CEDAR_2D_STENCIL_OP_BASE_H
#define CEDAR_2D_STENCIL_OP_BASE_H

#include <cedar/2d/grid_stencil.h>
#include <cedar/stencil_op_base.h>

namespace cedar { namespace cdr2 {

template <class grid_func, class registry>
class stencil_op_base : public ::cedar::stencil_op_base<grid_func, registry, grid_stencil>
{
public:
	stencil_op_base() {}
    stencil_op_base(len_t nx, len_t ny, bool intergrid=false) : gs(nx, ny, 1, intergrid) {}
	virtual real_t *data() { return gs.data(); }
	virtual grid_stencil & stencil() { return gs; }
	virtual const grid_stencil & stencil() const { return gs; }

protected:
	grid_stencil gs;
};

}}


#endif
