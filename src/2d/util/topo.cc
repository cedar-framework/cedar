#include <cmath>
#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"

#include <cedar/mpi/block_partition.h>
#include <cedar/decomp.h>

#include <cedar/2d/util/topo.h>

namespace cedar { namespace cdr2 { namespace util {

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


topo_ptr create_topo(MPI_Comm comm, int npx, int npy, len_t nx, len_t ny)
{
	int rank;

	auto igrd = std::make_shared<std::vector<len_t>>(NBMG_pIGRD);
	auto grid = std::make_shared<grid_topo>(igrd, 0, 1);

	grid->comm = comm;

	grid->nproc(0) = npx;
	grid->nproc(1) = npy;
	grid->nproc(2) = 1;

	MPI_Comm_rank(grid->comm, &rank);

	grid->coord(0) = rank % grid->nproc(0);
	grid->coord(1) = rank / grid->nproc(0);

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


topo_ptr create_topo(config & conf, MPI_Comm comm)
{
	auto islocal = conf.get<bool>("grid.local", true);
	auto ndofs = conf.getvec<len_t>("grid.n");
	auto nx = ndofs[0];
	auto ny = ndofs[1];
	topo_ptr grid;

	if (islocal) {
		auto nprocs = conf.getvec<int>("grid.np");
		int npx = 0;
		int npy = 0;
		if (nprocs.size() >= 2) {
			npx = nprocs[0];
			npy = nprocs[1];
		}
		if (npx == 0 or npy == 0) {
			grid = create_topo(comm, nx, ny);
		} else {
			int size;
			MPI_Comm_size(comm, &size);
			assert(size == npx*npy);
			grid = create_topo(comm, npx, npy, nx, ny);
		}
		log::status << "Running local solve" << std::endl;
	} else {
		grid = create_topo_global(comm, nx, ny);
		log::status << "Running global solve" << std::endl;
	}

	return grid;
}


topo_ptr model_topo(int nprocx, int nprocy, len_t nx, len_t ny)
{
	auto igrd = std::make_shared<std::vector<len_t>>(NBMG_pIGRD);
	auto grid = std::make_shared<grid_topo>(igrd, 0, 1);

	grid->nproc(0) = nprocx;
	grid->nproc(1) = nprocy;
	grid->nproc(2) = 1;

	grid->nlocal(0) = std::ceil(nx / static_cast<double>(grid->nproc(0)));
	grid->nlocal(1) = std::ceil(ny / static_cast<double>(grid->nproc(1)));

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


topo_ptr create_topo_global(MPI_Comm comm, len_t ngx, len_t ngy)
{
	int rank, size;

	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	auto igrd = std::make_shared<std::vector<len_t>>(NBMG_pIGRD);
	auto grid = std::make_shared<grid_topo>(igrd, 0, 1);

	grid->comm = comm;

	auto decomp = grid_decomp<2>(std::array<len_t,2>({{ngx,ngy}}), size);

	grid->nproc(0) = decomp[0];
	grid->nproc(1) = decomp[1];
	grid->nproc(2) = 1;

	assert(static_cast<unsigned int>(size) == grid->nproc());

	grid->coord(0) = rank % grid->nproc(0);
	grid->coord(1) = rank / grid->nproc(0);

	auto xpart = block_partition(ngx, grid->nproc(0));
	auto ypart = block_partition(ngy, grid->nproc(1));

	grid->nglobal(0) = ngx + 2;
	grid->nglobal(1) = ngy + 2;

	grid->is(0) = xpart.low(grid->coord(0)) + 1;
	grid->is(1) = ypart.low(grid->coord(1)) + 1;

	grid->nlocal(0) = xpart.size(grid->coord(0)) + 2;
	grid->nlocal(1) = ypart.size(grid->coord(1)) + 2;

	grid->dimxfine.resize(grid->nproc(0));
	grid->dimyfine.resize(grid->nproc(1));
	for (auto i : range(grid->nproc(0))) {
		grid->dimxfine[i] = xpart.size(i);
	}

	for (auto j : range(grid->nproc(1))) {
		grid->dimyfine[j] = ypart.size(j);
	}

	return grid;
}

}}}
