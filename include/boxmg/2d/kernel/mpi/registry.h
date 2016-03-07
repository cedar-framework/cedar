#ifndef BOXMG_2D_KERNEL_MPI_REGISTRY_H
#define BOXMG_2D_KERNEL_MPI_REGISTRY_H


#include "boxmg/types.h"
#include "boxmg/mpi_registry.h"
#include "boxmg/cycle/types.h"
#include "boxmg/kernel_name.h"

#include "boxmg/2d/relax_stencil.h"
#include "boxmg/2d/mpi/grid_func.h"
#include "boxmg/mpi/grid_topo.h"
#include "boxmg/2d/mpi/grid_func.h"



namespace boxmg { namespace bmg2d {
			class solver;
}}

namespace boxmg { namespace bmg2d { namespace mpi {
			class stencil_op;
			class redist_solver;
		}
		namespace inter {
			namespace mpi {
				class prolong_op;
			}
		}
	}
}


namespace boxmg { namespace bmg2d { namespace kernel { namespace mpi {
namespace mpi = boxmg::bmg2d::mpi;
class registry : public mpi_registry<mpi::stencil_op, relax_stencil, inter::mpi::prolong_op, mpi::grid_func, mpi::redist_solver, solver>
{
};

}}}}


#endif
