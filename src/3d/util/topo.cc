#include "boxmg/2d/ftn/mpi/BMG_workspace_c.h"

#include "boxmg/3d/util/topo.h"

namespace boxmg { namespace bmg3 { namespace util {

topo_ptr create_topo(MPI_Comm comm, len_t nx, len_t ny, len_t nz)
{
	int rank, size;
	int periodic[3];

	auto igrd = std::make_shared<std::vector<len_t>>(NBMG_pIGRD);

	auto grid = std::make_shared<grid_topo>(igrd, 0, 1);

	periodic[0] = periodic[1] = periodic[2] = 0;

	MPI_Comm_size(comm, &size);
	grid->nproc(0) = grid->nproc(1) = grid->nproc(2) = 0;
	MPI_Dims_create(size, 3, &grid->nproc(0));
	MPI_Cart_create(comm, 3, &grid->nproc(0), periodic, 1, &grid->comm);
	MPI_Comm_rank(grid->comm, &rank);
	MPI_Cart_coords(grid->comm, rank, 3, &grid->coord(0));

	auto tmp = grid->coord(0);
	grid->coord(0) = grid->coord(2);
	grid->coord(2) = tmp;

	tmp = grid->nproc(0);
	grid->nproc(0) = grid->nproc(2);
	grid->nproc(2) = tmp;

	grid->is(0) = grid->coord(0) * nx + 1;
	grid->nlocal(0) = nx;
	grid->is(1) = grid->coord(1) * ny + 1;
	grid->nlocal(1) = ny;
	grid->is(2) = grid->coord(2) * nz + 1;
	grid->nlocal(2) = nz;

	grid->nglobal(0) = nx*grid->nproc(0) + 2;
	grid->nglobal(1) = ny*grid->nproc(1) + 2;
	grid->nglobal(2) = nz*grid->nproc(2) + 2;

	grid->nlocal(0) += 2;
	grid->nlocal(1) += 2;
	grid->nlocal(2) += 2;

	// {
	// 	int rank;
	// 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// 	printf("[%d] %d %d %d -> %u %u %u : %u %u %u ==== %u %u %u\n", rank+1,
	// 	       grid->coord(0)+1, grid->coord(1)+1, grid->coord(2)+1,
	// 	       grid->nlocal(0), grid->nlocal(1), grid->nlocal(2),
	// 	       grid->nglobal(0), grid->nglobal(1), grid->nlocal(2),
	// 	       grid->is(0), grid->is(1), grid->is(2));
	// }

	return grid;
}

}}}
