#include <iostream>
#include <memory>
#include <math.h>

#include <cedar/types.h>
#include <cedar/2d/grid_func.h>
#include <cedar/2d/stencil_op.h>
#include <cedar/2d/solver.h>
#include <cedar/2d/gallery.h>

using namespace cedar;
using namespace cedar::cdr2;

static void set_problem(grid_func & b)
{
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


int main(int argc, char *argv[])
{
	config conf;
	cedar::init(conf);
	auto ndofs = conf.getvec<len_t>("grid.n");
	auto nx = ndofs[0];
	auto ny = ndofs[1];

	auto so = gallery::fe(nx, ny);
	grid_func b(nx, ny);
	set_problem(b);

	solver<nine_pt> bmg(so);
	auto x = grid_func::zeros(nx, ny);

	auto offload = conf.get<bool>("solver.offload");
	if (offload) {
		memory::prefetch(x.data(), x.size());
		memory::prefetch(b.data(), b.size());
		memory::sync();
	}
	bmg.solve(b, x);
	auto openmp = conf.get<bool>("solver.openmp");
	if (offload)
		timer_save("offload-" + std::to_string(nx*ny) + ".json");
	else if (openmp)
		timer_save("openmp-" + std::to_string(nx*ny) + ".json");
	else
		timer_save("serial-" + std::to_string(nx*ny) + ".json");

	log::status << "Finished Test" << std::endl;

	return 0;
}
