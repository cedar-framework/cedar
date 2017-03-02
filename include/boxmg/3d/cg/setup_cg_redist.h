#ifndef BOXMG_3D_KERNEL_SETUP_CG_REDIST_H
#define BOXMG_3D_KERNEL_SETUP_CG_REDIST_H

#include <memory>

#include <boxmg/3d/mpi/solver.h>
#include <boxmg/3d/mpi/redist_solver.h>
#include <boxmg/3d/mpi/stencil_op.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls {
	namespace mpi = boxmg::bmg3::mpi;
	void setup_cg_redist(const mpi::stencil_op & so,
	                     std::shared_ptr<config::reader> conf,
	                     std::shared_ptr<mpi::redist_solver> * slv,
	                     std::vector<int> & nblocks);
}

}}}
#endif
