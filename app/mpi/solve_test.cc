#include <mpi.h>
#include <iostream>
#include <memory>
#include <math.h>

#include <boxmg/types.h>
#include <boxmg/2d/mpi/grid_func.h>
#include <boxmg/2d/util/topo.h>
#include <boxmg/2d/util/mpi_grid.h>
#include <boxmg/2d/mpi/solver.h>


int main(int argc, char *argv[])
{
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	config::reader conf;

	auto islocal = conf.get<bool>("grid.local", true);
	auto nx = conf.get<len_t>("grid.nx", 9);
	auto ny = conf.get<len_t>("grid.ny", 9);
	topo_ptr grid;
	if (islocal) {
		grid = bmg2d::util::create_topo(MPI_COMM_WORLD, nx, ny);
		log::status << "Running local solve" << std::endl;
	} else {
		grid = bmg2d::util::create_topo_global(MPI_COMM_WORLD, nx, ny);
		log::status << "Running global solve" << std::endl;
	}

	auto so = mpi::stencil_op(grid);

	const double pi = M_PI;

	real_t hx = (1.0/(so.grid().nglobal(0)-1));
	real_t hy = (1.0/(so.grid().nglobal(1)-1));
	real_t h2 = hx * hy;

	grid_stencil & sten = so.stencil();
	sten.five_pt() = true;

	mpi::grid_func b(so.grid_ptr());

	int rank;
	MPI_Comm_rank(so.grid().comm, &rank);

	real_t y = so.grid().is(1)*hy;
	for (auto j: sten.range(1)) {
		real_t x = so.grid().is(0)*hx;
		for (auto i: sten.range(0)) {
			sten(i,j,dir::E) = 1;
			sten(i,j,dir::N) = 1;
			sten(i,j,dir::C) = 4;
			sten(i,j,dir::S) = 1;
			sten(i,j,dir::W) = 1;

			b(i,j) = 8*(pi*pi)*sin(2*pi*x)*sin(2*pi*y) * h2;
			x += hx;
		}
		y += hy;
	}

	if (util::mpi::has_boundary(so.grid(), dir::N))
		for (auto idx: sten.boarder(dir::N)) sten(idx, dir::N) = 0;
	if (util::mpi::has_boundary(so.grid(), dir::S))
		for (auto idx: sten.boarder(dir::S)) sten(idx, dir::S) = 0;
	if (util::mpi::has_boundary(so.grid(), dir::E))
		for (auto idx: sten.boarder(dir::E)) sten(idx, dir::E) = 0;
	if (util::mpi::has_boundary(so.grid(), dir::W))
	    for (auto idx: sten.boarder(dir::W)) sten(idx, dir::W) = 0;

	mpi::solver bmg(std::move(so));

	auto sol = bmg.solve(b);

	mpi::grid_func exact_sol(sol.grid_ptr());

	y = sol.grid().is(1)*hy;
	for (auto j: exact_sol.range(1)) {
		real_t x = sol.grid().is(0)*hx;
		for (auto i: exact_sol.range(0)) {
			exact_sol(i,j) = sin(2*pi*x) * sin(2*pi*y);
			x += hx;
		}
		y+= hy;
	}

	mpi::grid_func diff = exact_sol - sol;

	log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
