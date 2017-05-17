#ifndef CEDAR_2D_KERNEL_MPI_REGISTRY_H
#define CEDAR_2D_KERNEL_MPI_REGISTRY_H


#include "cedar/types.h"
#include "cedar/mpi_registry.h"
#include "cedar/cycle/types.h"
#include "cedar/kernel_name.h"

#include "cedar/2d/relax_stencil.h"
#include "cedar/2d/mpi/grid_func.h"
#include "cedar/mpi/grid_topo.h"
#include "cedar/2d/mpi/grid_func.h"



namespace cedar { namespace cdr2 {
			class solver;
}}

namespace cedar { namespace cdr2 { namespace mpi {
			class stencil_op;
			class redist_solver;
		}
		namespace inter {
			namespace mpi {
				class prolong_op;
				class restrict_op;
			}
		}
	}
}


namespace cedar { namespace cdr2 { namespace kernel { namespace mpi {
namespace mpi = cedar::cdr2::mpi;
class registry : public mpi_registry<mpi::stencil_op, relax_stencil, inter::mpi::prolong_op, inter::mpi::restrict_op, mpi::grid_func, mpi::redist_solver, solver>
{
};

}}}}


#endif
