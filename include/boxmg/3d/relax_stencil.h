#ifndef BOXMG_3D_RELAX_STENCIL_H
#define BOXMG_3D_RELAX_STENCIL_H

#include <boxmg/3d/types.h>
#include <boxmg/array.h>
#include <boxmg/grid_quantity.h>

namespace boxmg { namespace bmg3 {


class relax_stencil : public array<len_t, real_t,4>, public grid_quantity<len_t, 3>
{
public:
	relax_stencil() {};
	relax_stencil(len_t nx, len_t ny, len_t nz, unsigned int nghosts=1);

	using array<len_t,real_t,4>::index;
	using array<len_t,real_t,4>::operator();
};


}}

#endif
