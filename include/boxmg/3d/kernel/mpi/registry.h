#ifndef BOXMG_3D_KERNEL_MPI_REGISTRY_H
#define BOXMG_3D_KERNEL_MPI_REGISTRY_H

#include "boxmg/types.h"
#include "boxmg/mpi_registry.h"
#include "boxmg/cycle/types.h"
#include "boxmg/kernel_name.h"

#include "boxmg/3d/relax_stencil.h"
#include "boxmg/3d/mpi/grid_func.h"
#include "boxmg/mpi/grid_topo.h"
#include "boxmg/3d/mpi/grid_func.h"


namespace boxmg { namespace bmg3 {
			class solver;
}}

namespace boxmg { namespace bmg3 { namespace mpi {
			class solver;
			class stencil_op;
		}
		namespace inter {
			namespace mpi {
				class prolong_op;
			}
		}
	}
}


namespace boxmg { namespace bmg3 { namespace kernel { namespace mpi {
namespace mpi = boxmg::bmg3::mpi;
class registry : public mpi_registry<mpi::stencil_op, relax_stencil, inter::mpi::prolong_op, mpi::grid_func, mpi::solver, solver>
{
};

}}}}

#endif
