#include <memory>
#include <math.h>
#include <gtest/gtest.h>

#include <boxmg/types.h>
#include <boxmg/util/grid.h>
#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/solver.h>
#include <boxmg/3d/gallery.h>


static void set_problem(boxmg::bmg3::grid_func & b)
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

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


static void set_solution(boxmg::bmg3::grid_func & q)
{
	using namespace boxmg;

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


TEST(SerialPoisson3, Isotropic) {

	using namespace boxmg;
	using namespace boxmg::bmg3;

	auto nx = 200;
	auto ny = nx;
	auto nz = nx;

	auto so = gallery::poisson(nx, ny, nz);
	grid_func b(nx, ny, nz);

	set_problem(b);

	auto conf = std::make_shared<config::reader>("");
	solver bmg(std::move(so), conf);

	auto sol = bmg.solve(b);

	ASSERT_LT(std::abs(bmg.level(-1).res.lp_norm<2>()),
	          1e-8);

	grid_func exact_sol(sol.shape(0), sol.shape(1), sol.shape(2));

	set_solution(exact_sol);

	auto diff = exact_sol - sol;

	ASSERT_LT(std::abs(diff.inf_norm()),
	          1e-4);
}
