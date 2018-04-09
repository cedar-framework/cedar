#ifndef CEDAR_2D_TYPES_H
#define CEDAR_2D_TYPES_H

#include <cedar/solver_types.h>
#include <cedar/2d/base_types.h>
#include <cedar/2d/stencil_op.h>
#include <cedar/2d/grid_func.h>
#include <cedar/2d/restrict_op.h>
#include <cedar/2d/prolong_op.h>
#include <cedar/2d/relax_stencil.h>

namespace cedar
{
	namespace cdr2
	{
		using stypes = solver_types<
			stencil_op,
			five_pt,
			nine_pt,
			grid_func,
			prolong_op,
			restrict_op,
			relax_stencil>;
	}
}

#endif
