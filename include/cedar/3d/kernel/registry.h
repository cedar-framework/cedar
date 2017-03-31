#ifndef CEDAR_3D_KERNEL_REGISTRY_H
#define CEDAR_3D_KERNEL_REGISTRY_H


#include "cedar/kernel_registry.h"
#include "cedar/types.h"
#include "cedar/cycle/types.h"

#include "cedar/3d/relax_stencil.h"
#include "cedar/3d/grid_func.h"

namespace cedar { namespace cdr3 {
		class stencil_op;
		namespace inter {
			class prolong_op;
		}
	}
}

namespace cedar { namespace cdr3 { namespace kernel {

class registry : public kernel_registry<stencil_op, relax_stencil, inter::prolong_op, grid_func>
{
};

}}}


#endif
