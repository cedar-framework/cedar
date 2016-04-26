#include <mpi.h>
#include <iostream>

#include <boxmg/types.h>
#include <boxmg/3d/mpi/grid_func.h>
#include <boxmg/3d/util/topo.h>
#include <boxmg/3d/mpi/solver.h>

extern "C" {
	using namespace boxmg;
	void putf(real_t *so, real_t *qf,
	          len_t nlx, len_t nly, len_t lnz,
	          len_t ngx, len_t ngy, len_t ngz,
	          len_t igs, len_t jgs, len_t kgs,
	          real_t hx, real_t hy, real_t hz);
}

static void set_problem(boxmg::bmg3::mpi::stencil_op & so,
                        boxmg::bmg3::mpi::grid_func & b)
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

	real_t hx = (1.0/(so.grid().nglobal(0)-1));
	real_t hy = (1.0/(so.grid().nglobal(1)-1));
	real_t hz = (1.0/(so.grid().nglobal(2)-1));

	auto & sten = so.stencil();
	sten.five_pt() = true;
	auto & topo = so.grid();

	putf(so.data(), b.data(),
	     topo.nlocal(0)-2, topo.nlocal(1)-2, topo.nlocal(2)-2,
	     topo.nglobal(0)-2, topo.nglobal(1)-2, topo.nglobal(2)-2,
	     topo.is(0), topo.is(1), topo.is(2),
	     hx, hy, hz);
}


int main(int argc, char *argv[])
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

	int provided, rank;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	config::reader conf;

	auto islocal = conf.get<bool>("grid.local", true);
	auto nx = conf.get<len_t>("grid.nx", 31);
	auto ny = conf.get<len_t>("grid.ny", 31);
	auto nz = conf.get<len_t>("grid.nz", 31);
	topo_ptr grid;
	if (islocal) {
		auto npx = conf.get<int>("grid.npx", 0);
		auto npy = conf.get<int>("grid.npy", 0);
		auto npz = conf.get<int>("grid.npz", 0);
		if (npx == 0 or npy == 0 or npz == 0) {
			grid = bmg3::util::create_topo(MPI_COMM_WORLD, nx, ny, nz);
		} else {
			int size;
			MPI_Comm_size(MPI_COMM_WORLD, &size);
			assert(size == npx*npy*npz);
			grid = bmg3::util::create_topo(MPI_COMM_WORLD, npx, npy, npz,
			                              nx, ny, nz);

		}
		log::status << "Running local solve" << std::endl;
	} else {
		grid = bmg3::util::create_topo_global(MPI_COMM_WORLD, nx, ny, nz);
		log::status << "Running global solve" << std::endl;
	}
	//auto grid = util::create_topo(MPI_COMM_WORLD, nx, ny, nz);

	auto so = mpi::stencil_op(grid);
	mpi::grid_func b(so.grid_ptr());

	set_problem(so, b);
	mpi::solver bmg(std::move(so));
	auto x = bmg.solve(b);

	// {
	// 	std::ofstream ofile;
	// 	auto & topo = bmg.level(-1).A.grid();
	// 	ofile.open("op-" + std::to_string(topo.coord(0)+1) + "." +
	// 	           std::to_string(topo.coord(1)+1) + "." + std::to_string(topo.coord(2)+1) + ".txt", std::ios::out | std::ios::trunc | std::ios::binary);
	// 	ofile << bmg.level(0).A;
	// 	ofile.close();
	// }

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
