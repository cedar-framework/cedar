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

	auto offload = conf.get<bool>("solver.offload");
	auto openmp = conf.get<bool>("solver.openmp");
	auto kname = conf.get<std::string>("kernel");

	auto x = grid_func::ones(nx, ny);
	auto & flvl = bmg.levels.get(0);
	auto & clvl = bmg.levels.get(1);
	if (offload) {
		if (kname == "residual") {
			memory::prefetch(x.data(), x.size());
			memory::prefetch(b.data(), b.size());
			memory::prefetch(flvl.A.data(), flvl.A.size());
			memory::prefetch(flvl.res.data(), flvl.res.size());
		} else if (kname == "relax") {
			memory::prefetch(x.data(), x.size());
			memory::prefetch(b.data(), b.size());
			memory::prefetch(flvl.A.data(), flvl.A.size());
			memory::prefetch(flvl.SOR[0].data(), flvl.SOR[0].size());
		} else if (kname == "restrict") {
			memory::prefetch(clvl.P.data(), clvl.P.size());
			memory::prefetch(flvl.res.data(), flvl.res.size());
			memory::prefetch(clvl.b.data(), clvl.b.size());
		} else if (kname == "interp") {
			memory::prefetch(clvl.P.data(), clvl.P.size());
			memory::prefetch(x.data(), x.size());
			memory::prefetch(clvl.x.data(), clvl.x.size());
			memory::prefetch(flvl.res.data(), flvl.res.size());
			memory::prefetch(flvl.A.data(), flvl.A.size());
		}
		memory::sync();
	}
	auto kman = bmg.get_kernels();
	constexpr int nruns = 10;
	timer_begin("gpukern");
	if (kname == "residual") {
		for (int i = 0; i < nruns; i++) {
			kman->run<residual>(flvl.A, x, b, flvl.res);
		}
	} else if (kname == "relax") {
		for (int i = 0; i < nruns; i++) {
			kman->run<point_relax>(flvl.A, x, b, flvl.SOR[0], cycle::Dir::DOWN);
		}
	} else if (kname == "restrict") {
		for (int i = 0; i < nruns; i++) {
			kman->run<restriction>(clvl.R, flvl.res, clvl.b);
		}
	} else if (kname == "interp") {
		for (int i = 0; i < nruns; i++) {
			kman->run<interp_add>(clvl.P, clvl.x, flvl.res, x);
		}
	}
	timer_end("gpukern");

	if (offload)
		timer_save("offload-" + kname + "-" + std::to_string(nx*ny) + ".json");
	else if (openmp)
		timer_save("openmp-" + kname + "-" + std::to_string(nx*ny) + ".json");
	else
		timer_save("serial-" + kname + "-" + std::to_string(nx*ny) + ".json");

	log::status << "Finished Test" << std::endl;

	return 0;
}
