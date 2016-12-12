#include <memory>
#include <math.h>
#include <gtest/gtest.h>

#include <boxmg/types.h>
#include <boxmg/util/grid.h>
#include <boxmg/2d/grid_func.h>
#include <boxmg/2d/stencil_op.h>
#include <boxmg/2d/gallery.h>
#include <boxmg/2d/solver.h>


static void set_problem(boxmg::bmg2d::grid_func & b)
{
	using namespace boxmg;
	using namespace boxmg::bmg2d;

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


TEST(poisson, serial_2d) {
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	auto nx = 200;
	auto ny = nx;

	auto so = create_poisson(nx, ny);
	grid_func b(nx, ny);

	set_problem(b);

	config::reader conf("");
	solver bmg(std::move(so), std::move(conf));

	auto sol = bmg.solve(b);

	ASSERT_LT(std::abs(bmg.level(-1).res.lp_norm<2>()),
	          1e-8);

	grid_func exact_sol(sol.shape(0), sol.shape(1));

	set_solution(exact_sol);

	auto diff = exact_sol - sol;

	ASSERT_LT(std::abs(diff.inf_norm()),
	          1e-4);
}


TEST(poisson_x, serial_2d) {
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	auto nx = 800;
	auto ny = 200;

	auto so = create_poisson(nx, ny);
	grid_func b(nx, ny);

	set_problem(b);

	config::reader conf("");
	conf.set("solver.relaxation", "line-x");
	solver bmg(std::move(so), std::move(conf));

	auto sol = bmg.solve(b);

	ASSERT_LT(std::abs(bmg.level(-1).res.lp_norm<2>()),
	          1e-8);

	grid_func exact_sol(sol.shape(0), sol.shape(1));

	set_solution(exact_sol);

	auto diff = exact_sol - sol;

	ASSERT_LT(std::abs(diff.inf_norm()),
	          1e-4);
}


TEST(poisson_y, serial_2d) {
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	auto nx = 200;
	auto ny = 800;

	auto so = create_poisson(nx, ny);
	grid_func b(nx, ny);

	set_problem(b);

	config::reader conf("");
	conf.set("solver.relaxation", "line-y");
	solver bmg(std::move(so), std::move(conf));

	auto sol = bmg.solve(b);

	ASSERT_LT(std::abs(bmg.level(-1).res.lp_norm<2>()),
	          1e-8);

	grid_func exact_sol(sol.shape(0), sol.shape(1));

	set_solution(exact_sol);

	auto diff = exact_sol - sol;

	ASSERT_LT(std::abs(diff.inf_norm()),
	          1e-4);
}
