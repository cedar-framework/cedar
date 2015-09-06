#include <mpi.h>
#include <iostream>
#include <memory>
#include <math.h>

#include <boxmg-common.h>
#include <boxmg-2d.h>


int main(int argc, char *argv[])
{
	using namespace boxmg;
	using namespace boxmg::bmg2d;
	using namespace boxmg::bmg2d::core;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, NULL);

	auto nx = config::get<len_t>("grid.nx", 9);
	auto ny = config::get<len_t>("grid.ny", 9);
	auto grid = bmg2d::util::create_topo(MPI_COMM_WORLD, nx, ny);
	auto so = mpi::StencilOp(grid);

	const double pi = M_PI;

	real_t hx = (1.0/(so.grid().nglobal(0)-3));
	real_t hy = (1.0/(so.grid().nglobal(1)-3));
	real_t h2 = hx * hy;

	GridStencil & sten = so.stencil();
	sten.five_pt() = true;

	mpi::GridFunc b(so.grid_ptr());

	int rank;
	MPI_Comm_rank(so.grid().comm, &rank);

	real_t y = (so.grid().is(1)-1)*hy;
	for (auto j: sten.range(1)) {
		real_t x = (so.grid().is(0)-1)*hx;
		for (auto i: sten.range(0)) {
			log::debug << i << " " << j << std::endl;
			sten(i,j,Dir::E) = 1;
			sten(i,j,Dir::N) = 1;
			sten(i,j,Dir::C) = 4;
			sten(i,j,Dir::S) = 1;
			sten(i,j,Dir::W) = 1;

			b(i,j) = 8*(pi*pi)*sin(2*pi*x)*sin(2*pi*y) * h2;
			x += hx;
		}
		y += hy;
	}

	if (util::mpi::has_boundary(so.grid(), Dir::N)) {
		for (auto idx: sten.boarder(Dir::N)){
			sten(idx, Dir::C) = 1;
			sten(idx, Dir::S) = 0;
			sten(idx, Dir::E) = 0;
			sten(idx, Dir::W) = 0;
			sten(idx, Dir::N) = 0;
		}
	}
	if (util::mpi::has_boundary(so.grid(), Dir::S)) {
		for (auto idx: sten.boarder(Dir::S)){
			sten(idx, Dir::C) = 1;
			sten(idx, Dir::S) = 0;
			sten(idx, Dir::E) = 0;
			sten(idx, Dir::W) = 0;
			sten(idx, Dir::N) = 0;
		}
	}
	if (util::mpi::has_boundary(so.grid(), Dir::E)) {
		for (auto idx: sten.boarder(Dir::E)) {
			sten(idx, Dir::C) = 1;
			sten(idx, Dir::S) = 0;
			sten(idx, Dir::E) = 0;
			sten(idx, Dir::W) = 0;
			sten(idx, Dir::N) = 0;
		}
	}
	if (util::mpi::has_boundary(so.grid(), Dir::W)) {
		for (auto idx: sten.boarder(Dir::W)) {
			sten(idx, Dir::C) = 1;
			sten(idx, Dir::S) = 0;
			sten(idx, Dir::E) = 0;
			sten(idx, Dir::W) = 0;
			sten(idx, Dir::N) = 0;
		}
	}

	// std::ofstream ofile;
	// ofile.open("rhs-" + std::to_string(rank) + ".txt", std::ios::out | std::ios::trunc | std::ios::binary);
	// ofile << b;
	// ofile.close();

	solver::mpi::BoxMG bmg(std::move(so));

	auto sol = bmg.solve(b);

	auto &op = bmg.level(-2).A;
	auto & gt = op.grid();
	std::ofstream ofile;
	ofile.open("op" + std::to_string(gt.coord(0)) + "-" + std::to_string(gt.coord(1)) + ".txt", std::ios::out | std::ios::trunc | std::ios::binary);
	ofile << op;
	ofile.close();

	core::mpi::GridFunc exact_sol(sol.grid_ptr());

	y = (sol.grid().is(1)-1)*hy;
	for (auto j: exact_sol.range(1)) {
		real_t x = (sol.grid().is(0)-1)*hx;
		for (auto i: exact_sol.range(0)) {
			exact_sol(i,j) = sin(2*pi*x) * sin(2*pi*y);
			x += hx;
		}
		y+= hy;
	}

	auto diff = exact_sol - sol;

	log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
