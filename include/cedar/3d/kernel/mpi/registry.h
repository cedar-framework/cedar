#ifndef CEDAR_3D_KERNEL_MPI_REGISTRY_H
#define CEDAR_3D_KERNEL_MPI_REGISTRY_H

#include "cedar/types.h"
#include "cedar/mpi_registry.h"
#include "cedar/cycle/types.h"
#include "cedar/kernel_name.h"

#include "cedar/3d/relax_stencil.h"
#include "cedar/3d/mpi/grid_func.h"
#include "cedar/mpi/grid_topo.h"
#include "cedar/3d/mpi/grid_func.h"


namespace cedar { namespace cdr3 {
			class solver;
}}

namespace cedar { namespace cdr3 { namespace mpi {
			class solver;
			class redist_solver;
			class stencil_op;
		}
		namespace inter {
			namespace mpi {
				class prolong_op;
			}
		}
	}
}


namespace cedar { namespace cdr3 { namespace kernel { namespace mpi {
namespace mpi = cedar::cdr3::mpi;
class registry : public mpi_registry<mpi::stencil_op, relax_stencil, inter::mpi::prolong_op, mpi::grid_func, mpi::redist_solver, solver>
{
};

}}}}

#endif
