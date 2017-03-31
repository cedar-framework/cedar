#ifndef CEDAR_2D_KERNEL_REGISTRY_H
#define CEDAR_2D_KERNEL_REGISTRY_H


#include "cedar/kernel_registry.h"
#include "cedar/types.h"
#include "cedar/cycle/types.h"

#include "cedar/2d/relax_stencil.h"
#include "cedar/2d/grid_func.h"

namespace cedar { namespace cdr2 {
		class stencil_op;
		namespace inter {
			class prolong_op;
		}
	}
}

namespace cedar { namespace cdr2 { namespace kernel {

class registry : public kernel_registry<stencil_op, relax_stencil, inter::prolong_op, grid_func>
{
};

}}}


#endif
