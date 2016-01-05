#ifndef BOXMG_3D_STENCIL_OP_BASE_H
#define BOXMG_3D_STENCIL_OP_BASE_H

#include <boxmg/stencil_op_base.h>
#include <boxmg/3d/grid_stencil.h>

namespace boxmg { namespace bmg3 {

template <class grid_func, class registry>
class stencil_op_base : public ::boxmg::stencil_op_base<grid_func, registry, grid_stencil>
{
public:
	stencil_op_base() {}
    stencil_op_base(len_t nx, len_t ny, len_t nz, bool intergrid=false) : gs(nx, ny, nz, 1, intergrid) {}
	virtual real_t *data() { return gs.data(); }
	virtual grid_stencil & stencil() { return gs; }
	virtual const grid_stencil & stencil() const { return gs; }

protected:
	grid_stencil gs;

};

}}
#endif
