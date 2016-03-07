#ifndef BOXMG_2D_MPI_REDIST_SOLVER_H
#define BOXMG_2D_MPI_REDIST_SOLVER_H

#include <mpi.h>
#include <array>

#include <boxmg/2d/mpi/solver.h>
#include <boxmg/2d/stencil_op.h>
#include <boxmg/2d/mpi/halo.h>

namespace boxmg { namespace bmg2d { namespace mpi {

class redist_solver
{
public:
	using msg_ctx = boxmg::bmg2d::kernel::impls::MsgCtx;
	redist_solver(const stencil_op & so, std::array<int, 2> nblock);
	std::shared_ptr<grid_topo> redist_topo(const grid_topo & fine_topo, msg_ctx & ctx);

protected:
	MPI_Comm collcomm;
	std::array<int,2> nblock;
};

}}}

#endif
