#include "boxmg/2d/ftn/mpi/BMG_workspace_c.h"

#include <boxmg/mpi/block_partition.h>
#include <boxmg/2d/mpi/redist_solver.h>


using namespace boxmg;
using namespace boxmg::bmg2d::mpi;

redist_solver::redist_solver(const stencil_op & so, std::array<int, 2> nblock) : nblock(nblock)
{
	// Split communicator into collective processor blocks
	auto & topo = so.grid();
	msg_ctx * ctx = (msg_ctx*) so.halo_ctx;
	auto ctopo = redist_topo(topo, *ctx);
}


std::shared_ptr<grid_topo> redist_solver::redist_topo(const grid_topo & fine_topo, msg_ctx & ctx)
{
	// std::cout << fine_topo.coord(0) << " " << fine_topo.coord(1) << " => ("
	//           << fine_topo.nglobal(0) << ", " << fine_topo.nglobal(1) << ") ("
	//           << fine_topo.nlocal(0) << ", " << fine_topo.nlocal(1) << ")" << std::endl;
	auto igrd = std::make_shared<std::vector<len_t>>(NBMG_pIGRD);
	auto grid = std::make_shared<grid_topo>(igrd, 0, 1);

	grid->nproc(0) = nblock[0];
	grid->nproc(1) = nblock[1];

	grid->nglobal(0) = fine_topo.nglobal(0);
	grid->nglobal(1) = fine_topo.nglobal(1);

	grid->nlocal(0) = 0;
	grid->nlocal(1) = 0;

	grid->is(0) = 0;
	grid->is(1) = 0;

	// block mapping
	block_partition parti(fine_topo.nproc(0), nblock[0]);
	block_partition partj(fine_topo.nproc(1), nblock[1]);

	grid->coord(0) = parti.owner(fine_topo.coord(0));
	grid->coord(1) = partj.owner(fine_topo.coord(1));

	auto lowi = parti.low(grid->coord(0));
	auto highi = parti.high(grid->coord(0));
	auto lowj = partj.low(grid->coord(1));
	auto highj = partj.high(grid->coord(1));
	for (auto i = lowi; i <= highi; i++) {
		grid->nlocal(0) += ctx.cg_nlocal(0, ctx.proc_grid(i, lowj)) - 2; // remove ghosts
	}
	for (auto j = lowj; j <= highj; j++) {
		grid->nlocal(1) += ctx.cg_nlocal(1, ctx.proc_grid(lowi, j)) - 2; // remove ghosts
	}
	for (auto i = 0; i < lowi; i++) {
		grid->is(0) += ctx.cg_nlocal(0, ctx.proc_grid(i, lowj)) - 2; // remove ghosts
	}
	grid->is(0)++; // 1 based indexing
	for (auto j = 0; j < lowj; j++) {
		grid->is(1) += ctx.cg_nlocal(0, ctx.proc_grid(lowi, j)) - 2; // remove ghosts
	}
	grid->is(1)++;

	// add ghosts
	grid->nlocal(0) += 2;
	grid->nlocal(1) += 2;

	int color = grid->coord(0) + grid->nproc(0)*grid->coord(1);
	int key = (fine_topo.coord(0) - lowi) + (fine_topo.coord(1) - lowj)*parti.size(grid->coord(0));
	MPI_Comm_split(fine_topo.comm, color, key, &this->collcomm);

	MPI_Comm_split(fine_topo.comm, key, color, &grid->comm);

	// if (fine_topo.coord(0) == 0 and fine_topo.coord(1) == 0) {
	// 	std::cout << lowi << " " << highi << std::endl;
	// 	std::cout << lowj << " " << highj << std::endl;
	// 	std::cout << ctx.cg_nlocal(0, ctx.proc_grid(1,0)) - 2 << std::endl;
	// }

	// std::cout << grid->coord(0) << " " << grid->coord(1) << " => "
	//           << grid->nlocal(0)-2 << " "  << grid->nlocal(1)-2 << std::endl;

	std::cout << fine_topo.coord(0) << " " << fine_topo.coord(1) << " => ("
	          << color << " " << key << ")" << std::endl;

	MPI_Barrier(fine_topo.comm);
	MPI_Abort(fine_topo.comm, 0);

	return grid;
}
