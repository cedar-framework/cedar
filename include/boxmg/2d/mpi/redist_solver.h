#ifndef BOXMG_2D_MPI_REDIST_SOLVER_H
#define BOXMG_2D_MPI_REDIST_SOLVER_H

#include <mpi.h>
#include <array>

#include <boxmg/config/reader.h>
#include <boxmg/mpi/redist_comms.h>
#include <boxmg/2d/mpi/solver.h>
#include <boxmg/2d/stencil_op.h>
#include <boxmg/2d/mpi/halo.h>

namespace boxmg { namespace bmg2d { namespace mpi {

class redist_solver
{
public:
	using msg_ctx = boxmg::bmg2d::kernel::impls::MsgCtx;
	redist_solver(const stencil_op & so, std::shared_ptr<config::reader> conf, std::array<int, 2> nblock);
	void solve(const grid_func & b, grid_func & x);

protected:
	bool redundant;
	std::unique_ptr<solver> slv;
	redist_comms rcomms;
	int block_id; /** id within a block */
	int block_num; /** which block */
	std::array<int,2> nblock;
	bool active;
	int recv_id;
	std::vector<int> send_ids;
	array<len_t, len_t, 1> nbx; /** # d.o.f. for each processor in my block */
	array<len_t, len_t, 1> nby; /** # d.o.f. for each processor in my block */
	MPI_Fint msg_comm;
	grid_func b_redist;
	grid_func x_redist;

	std::shared_ptr<grid_topo> redist_topo(const grid_topo & fine_topo, msg_ctx & ctx);
	stencil_op redist_operator(const stencil_op & so, topo_ptr topo);
	void gather_rhs(const grid_func & b);
	void scatter_sol(grid_func & x);
};

}}}

#endif
