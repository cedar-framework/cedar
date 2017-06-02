#ifndef CEDAR_2D_KERNEL_SETUP_CG_REDIST_H
#define CEDAR_2D_KERNEL_SETUP_CG_REDIST_H

#include <memory>

#include <cedar/kernel_params.h>
#include <cedar/2d/mpi/stencil_op.h>

namespace cedar { namespace cdr2 { namespace mpi {
			class redist_solver;
		}
	}
}

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls {
	namespace mpi = cedar::cdr2::mpi;
	template <class sten>
	void setup_cg_redist(const kernel_params & params,
	                     const mpi::stencil_op<sten> & so,
	                     std::shared_ptr<config::reader> conf,
	                     std::shared_ptr<mpi::redist_solver> * slv,
	                     std::vector<int> & nblocks);
}

}}}
#endif

