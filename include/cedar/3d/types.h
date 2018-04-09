#ifndef CEDAR_3D_TYPES_H
#define CEDAR_3D_TYPES_H

#include <cedar/3d/base_types.h>
#include <cedar/solver_types.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/grid_func.h>
#include <cedar/3d/restrict_op.h>
#include <cedar/3d/prolong_op.h>
#include <cedar/3d/relax_stencil.h>

namespace cedar { namespace cdr3 {
		using stypes = solver_types<
			stencil_op,
			seven_pt,
			xxvii_pt,
			grid_func,
			prolong_op,
			restrict_op,
			relax_stencil>;
}}

#endif
