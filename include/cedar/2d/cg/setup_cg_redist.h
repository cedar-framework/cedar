#ifndef CEDAR_2D_KERNEL_SETUP_CG_REDIST_H
#define CEDAR_2D_KERNEL_SETUP_CG_REDIST_H

#include <memory>

#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/mpi/redist_solver.h>
#include <cedar/2d/mpi/stencil_op.h>

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls {
	namespace mpi = cedar::cdr2::mpi;
	void setup_cg_redist(const mpi::stencil_op & so,
	                     std::shared_ptr<config::reader> conf,
	                     std::shared_ptr<mpi::redist_solver> * slv,
	                     std::vector<int> & nblocks);
}

}}}
#endif

