#include <cmath>
#include "boxmg/2d/ftn/mpi/BMG_workspace_c.h"

#include "boxmg/2d/util/topo.h"

namespace boxmg { namespace bmg2d { namespace util {

topo_ptr create_topo(MPI_Comm comm, len_t nx, len_t ny)
{
	int rank, size;
	int periodic[3];

	auto igrd = std::make_shared<std::vector<len_t>>(NBMG_pIGRD);

	auto grid = std::make_shared<grid_topo>(igrd, 0, 1);

	periodic[0] = periodic[1] = periodic[2] = 0;

	MPI_Comm_size(comm, &size);
	grid->nproc(0) = grid->nproc(1) = 0;
	grid->nproc(2) = 1;
	MPI_Dims_create(size, 2, &grid->nproc(0));
	MPI_Cart_create(comm, 2, &grid->nproc(0), periodic, 1, &grid->comm);
	MPI_Comm_rank(grid->comm, &rank);
	MPI_Cart_coords(grid->comm, rank, 2, &grid->coord(0));

	auto tmp = grid->coord(0);
	grid->coord(0) = grid->coord(1);
	grid->coord(1) = tmp;

	tmp = grid->nproc(0);
	grid->nproc(0) = grid->nproc(1);
	grid->nproc(1) = tmp;

	grid->is(0) = grid->coord(0) * nx + 1;
	grid->nlocal(0) = nx;
	grid->is(1) = grid->coord(1) * ny + 1;
	grid->nlocal(1) = ny;

	grid->nglobal(0) = nx*grid->nproc(0) + 2;
	grid->nglobal(1) = ny*grid->nproc(1) + 2;

	grid->nlocal(0) += 2;
	grid->nlocal(1) += 2;

	// printf("%d %d -> %u %u : %u %u ==== %u %u\n", grid->coord(0), grid->coord(1),
	//        grid->nlocal(0), grid->nlocal(1),
	//        grid->nglobal(0), grid->nglobal(1),
	//        grid->is(0), grid->is(1));

	return grid;
}


topo_ptr create_topo(int np, len_t nx, len_t ny)
{
	auto igrd = std::make_shared<std::vector<len_t>>(NBMG_pIGRD);
	auto grid = std::make_shared<grid_topo>(igrd, 0, 1);

	grid->nproc(0) = grid->nproc(1) = 0;
	grid->nproc(2) = 1;

	MPI_Dims_create(np, 2, &grid->nproc(0));

	auto tmp = grid->nproc(0);
	grid->nproc(0) = grid->nproc(1);
	grid->nproc(1) = tmp;
	grid->nproc(2) = 1;

	grid->nlocal(0) = std::ceil(nx / grid->nproc(0));
	grid->nlocal(1) = std::ceil(ny / grid->nproc(1));

	grid->nglobal(0) = nx;
	grid->nglobal(1) = ny;

	return grid;
}


topo_ptr coarsen_topo(topo_ptr topof)
{
	auto igrd = std::make_shared<std::vector<len_t>>(NBMG_pIGRD);
	auto grid = std::make_shared<grid_topo>(igrd, 0, 1);

	grid->nproc(0) = topof->nproc(0);
	grid->nproc(1) = topof->nproc(1);
	grid->nproc(2) = topof->nproc(2);

	grid->nglobal(0) = (topof->nglobal(0) - 1) / 2 + 2;
	grid->nglobal(1) = (topof->nglobal(1) - 1) / 2 + 2;

	grid->nlocal(0) = (topof->nlocal(0) - 1) / 2 + 2;
	grid->nlocal(1) = (topof->nlocal(1) - 1) / 2 + 2;

	return grid;
}

}}}
