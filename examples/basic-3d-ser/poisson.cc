#include <math.h>
#include <memory>
#include <iostream>

#include <cedar/types.h>
#include <cedar/3d/gallery.h>
#include <cedar/3d/solver.h>

using namespace cedar;
using namespace cedar::cdr3;

static void set_problem(grid_func & b)
{
	const double pi = M_PI;

	auto rhs = [pi](real_t x, real_t y, real_t z) {
		return 12*(pi*pi)*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z);
	};

	b.set(0);

	real_t hx = 1.0 / (b.len(0) - 1);
	real_t hy = 1.0 / (b.len(1) - 1);
	real_t hz = 1.0 / (b.len(2) - 1);

	real_t h2 = hx*hy*hz;

	for (auto k : b.range(2)) {
		for (auto j : b.range(1)) {
			for (auto i : b.range(0)) {
				real_t x = i*hx;
				real_t y = j*hy;
				real_t z = k*hz;

				b(i,j,k) = rhs(x,y,z) * h2;
			}
		}
	}
}


static void set_solution(grid_func & q)
{
	const double pi = M_PI;

	auto sol = [pi](real_t x, real_t y, real_t z) {
		return sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z);
	};

	real_t hx = 1.0 / (q.len(0) - 1);
	real_t hy = 1.0 / (q.len(1) - 1);
	real_t hz = 1.0 / (q.len(2) - 1);

	for (auto k : q.range(2)) {
		for (auto j : q.range(1)) {
			for (auto i : q.range(0)) {
				real_t x = i*hx;
				real_t y = j*hy;
				real_t z = k*hz;

				q(i,j,k) = sol(x,y,z);
			}
		}
	}
}


int main(int argc, char *argv[])
{
	log::status << "Beginning test" << std::endl;

	config conf;

	auto ndofs = conf.getvec<len_t>("grid.n");
	auto nx = ndofs[0];
	auto ny = ndofs[1];
	auto nz = ndofs[2];

	auto so = gallery::poisson(nx, ny, nz);
	grid_func b(nx, ny, nz);
	set_problem(b);

	solver<seven_pt> bmg(so);

	auto sol = bmg.solve(b);

	grid_func exact_sol(nx, ny, nz);

	set_solution(exact_sol);

	auto diff = exact_sol - sol;

	log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	log::status << "Finished test" << std::endl;

	return 0;
}
