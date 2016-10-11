#include <iostream>
#include <memory>
#include <math.h>

#include <boxmg/types.h>
#include <boxmg/util/grid.h>
#include <boxmg/2d/grid_func.h>
#include <boxmg/2d/stencil_op.h>
#include <boxmg/2d/solver.h>


static void set_problem(boxmg::bmg2d::stencil_op & so, boxmg::bmg2d::grid_func & b)
{
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	const double pi = M_PI;

	auto rhs = [pi](real_t x, real_t y) {
		return 8*(pi*pi)*sin(2*pi*x)*sin(2*pi*y);
	};

	grid_stencil & sten = so.stencil();
	sten.five_pt() = true;

	b.set(0);
	sten.set(0);

	real_t hx = 1.0/(b.len(0)-1);
	real_t hy = 1.0/(b.len(1)-1);
	real_t h2 = hx*hy;
	real_t xh = hy/hx;
	real_t yh = hx/hy;
	len_t i1 = b.len(0);
	len_t j1 = b.len(1);

	for (auto j : range<len_t>(2, j1)) {
		for (auto i : range<len_t>(1, i1)) {
			sten(i,j,dir::S) = 1.0 * yh;
		}
	}

	for (auto j : range<len_t>(1, j1)) {
		for (auto i : range<len_t>(2, i1)) {
			sten(i,j,dir::W) = 1.0 * xh;
		}
	}

	for (auto j : b.range(1)) {
		for (auto i : b.range(0)) {
			real_t x = i*hx;
			real_t y = j*hy;

			sten(i,j,dir::C) = 2*xh + 2*yh;
			b(i,j) = rhs(x,y) * h2;
		}
	}
}


static void set_solution(boxmg::bmg2d::grid_func & q)
{
	using namespace boxmg;

	const double pi = M_PI;

	auto sol = [pi](real_t x, real_t y) {
		return sin(2*pi*x)*sin(2*pi*y);
	};

	real_t hx = 1.0/(q.len(0)-1);
	real_t hy = 1.0/(q.len(1)-1);
	for (auto j : q.grange(1)) {
		for (auto i : q.grange(0)) {
			real_t x = i*hx;
			real_t y = j*hy;

			q(i,j) = sol(x, y);
		}
	}
}


int main(int argc, char *argv[])
{
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	const double pi = M_PI;

	config::reader conf;
	auto ndofs = conf.getvec<len_t>("grid.n");
	auto nx = ndofs[0];
	auto ny = ndofs[1];

	real_t h2 = (1.0/(nx+1)) * (1.0/(ny+1));

	auto so = stencil_op(nx,ny);
	grid_func b(nx, ny);

	set_problem(so, b);

	solver bmg(std::move(so));

	auto sol = bmg.solve(b);

	grid_func exact_sol(sol.shape(0), sol.shape(1));

	std::ofstream dfile;
	dfile.open("sol.txt", std::ios::out | std::ios::trunc | std::ios::binary);

	set_solution(exact_sol);

	auto diff = exact_sol - sol;
	dfile << sol << std::endl;

	dfile.close();

	log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	log::status << "Finished Test" << std::endl;

	return 0;
}
