#ifndef BOXMG_3D_KERNEL_REGISTRY_H
#define BOXMG_3D_KERNEL_REGISTRY_H


#include "boxmg/kernel_registry.h"
#include "boxmg/types.h"
#include "boxmg/cycle/types.h"

#include "boxmg/3d/relax_stencil.h"
#include "boxmg/3d/grid_func.h"

namespace boxmg { namespace bmg3 {
		class stencil_op;
		namespace inter {
			class prolong_op;
		}
	}
}

namespace boxmg { namespace bmg3 { namespace kernel {

class registry : public kernel_registry<stencil_op, relax_stencil, inter::prolong_op, grid_func>
{
};

}}}


#endif
