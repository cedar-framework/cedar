#include <iostream>
#include <memory>
#include <random>
#include <math.h>

#include <cedar/types.h>
#include <cedar/2d/grid_func.h>
#include <cedar/2d/stencil_op.h>
#include <cedar/2d/solver.h>
#include <cedar/2d/gallery.h>

using namespace cedar;
using namespace cedar::cdr2;

static void set_random(grid_func & x)
{
	using namespace cedar;

	std::mt19937 gen;
	gen.seed(0);
	std::uniform_real_distribution<real_t> dis;

	for (int j = 1; j < x.len(1); j++) {
		for (int i = 1; i < x.len(0); i++) {
			x(i,j) = dis(gen);
		}
	}
}

template<class sten>
static void set_random(stencil_op<sten> & so)
{
	using namespace cedar;

	std::mt19937 gen;
	gen.seed(0);
	std::uniform_real_distribution<real_t> dis;

	so.set(0);

	for (int j = 1; j < so.len(1); j++) {
		for (int i = 1; i < so.len(0); i++) {
			for (auto k : range(stencil_ndirs<sten>::value)) {
				so(i,j,static_cast<sten>(k)) = dis(gen);
			}
		}
	}
}


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
	timer_save(offload ? "offload.json" : "system.json");

	log::status << "Finished Test" << std::endl;

	return 0;
}
