#include "relax_stencil.h"

using namespace boxmg::bmg2d;

relax_stencil::relax_stencil(len_t nx, len_t ny, unsigned int nghosts) : Array(nx,ny,nghosts,false)
{
	len_t size;
	nx += nghosts*2;
	ny += nghosts*2;

	size = nx*ny*2;

	data_.resize(size);
}
