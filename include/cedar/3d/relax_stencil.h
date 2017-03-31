#ifndef CEDAR_3D_RELAX_STENCIL_H
#define CEDAR_3D_RELAX_STENCIL_H

#include <cedar/3d/types.h>
#include <cedar/array.h>
#include <cedar/grid_quantity.h>

namespace cedar { namespace cdr3 {


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
