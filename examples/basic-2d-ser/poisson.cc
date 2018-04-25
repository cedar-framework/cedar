#include <iostream>
#include <memory>
#include <math.h>

#include <cedar/types.h>
#include <cedar/util/grid.h>
#include <cedar/2d/grid_func.h>
#include <cedar/2d/stencil_op.h>
#include <cedar/2d/solver.h>
#include <cedar/2d/gallery.h>



static void set_problem(cedar::cdr2::grid_func & b)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	const double pi = M_PI;

	auto rhs = [pi](real_t x, real_t y) {
		return 8*(pi*pi)*sin(2*pi*x)*sin(2*pi*y);
	};

	b.set(0);

	real_t hx = 1.0/(b.len(0)-1);
	real_t hy = 1.0/(b.len(1)-1);
	real_t h2 = hx*hy;

	for (auto j : b.range(1)) {
		for (auto i : b.range(0)) {
			real_t x = i*hx;
			real_t y = j*hy;

			b(i,j) = rhs(x,y) * h2;
		}
	}
}


static void set_solution(cedar::cdr2::grid_func & q)
{
	using namespace cedar;

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
	using namespace cedar;
	using namespace cedar::cdr2;

	config conf;
	auto ndofs = conf.getvec<len_t>("grid.n");
	auto nx = ndofs[0];
	auto ny = ndofs[1];

	auto so = gallery::poisson(nx, ny);
	grid_func b(nx, ny);

	set_problem(b);

	solver<five_pt> bmg(so);

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
