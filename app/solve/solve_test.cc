#include <iostream>
#include <memory>
#include <math.h>

#include <boxmg/types.h>
#include <boxmg/util/grid.h>
#include <boxmg/2d/grid_func.h>
#include <boxmg/2d/stencil_op.h>
#include <boxmg/2d/solver.h>


int main(int argc, char *argv[])
{
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	const double pi = M_PI;

	auto nx = config::get<len_t>("grid.nx", 301);
	auto ny = config::get<len_t>("grid.ny", 301);

	real_t h2 = (1.0/(nx+1)) * (1.0/(ny+1));

	auto so = stencil_op(nx,ny);
	grid_stencil & sten = so.stencil();
	sten.five_pt() = true;

	grid_func b(nx, ny);
	auto xs = linspace(0,1,nx+2);
	auto ys = linspace(0,1,ny+2);

	auto y = ys.begin();
	y++;
	for (auto j: sten.range(1)) {
		auto x = xs.begin();
		x++;
		for (auto i: sten.range(0)) {
			sten(i,j,dir::E) = 1;
			sten(i,j,dir::N) = 1;
			sten(i,j,dir::C) = 4;
			sten(i,j,dir::S) = 1;
			sten(i,j,dir::W) = 1;

			b(i,j) = 8*(pi*pi)*sin(2*pi*(*x))*sin(2*pi*(*y)) * h2;
			x++;
		}
		y++;
	}

	for (auto idx: sten.boarder(dir::N)) sten(idx, dir::N) = 0;
	for (auto idx: sten.boarder(dir::S)) sten(idx, dir::S) = 0;
	for (auto idx: sten.boarder(dir::E)) sten(idx, dir::E) = 0;
	for (auto idx: sten.boarder(dir::W)) sten(idx, dir::W) = 0;

	solver bmg(std::move(so));

	auto sol = bmg.solve(b);

	grid_func exact_sol(sol.shape(0), sol.shape(1));

	std::ofstream dfile;
	dfile.open("sol.txt", std::ios::out | std::ios::trunc | std::ios::binary);

	y = ys.begin();
	for (auto j: exact_sol.grange(1)) {
		auto x = xs.begin();
		for (auto i: exact_sol.grange(0)) {
			exact_sol(i,j) = sin(2*pi*(*x)) * sin(2*pi*(*y));
			x++;
		}
		y++;
	}


	auto diff = exact_sol - sol;
	dfile << sol << std::endl;

	dfile.close();

	log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	log::status << "Finished Test" << std::endl;

	return 0;
}
