#include <gtest/gtest.h>
#include <mpi.h>

#include <cedar/2d/mpi/gallery.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/mpi/solver.h>


static void set_problem(cedar::cdr2::mpi::grid_func & b, std::function<cedar::real_t(cedar::real_t,cedar::real_t)> rhs)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	auto & topo = b.grid();

	b.set(0);

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);

	real_t hx = 1.0 / (topo.nglobal(0) - 1);
	real_t hy = 1.0 / (topo.nglobal(1) - 1);

	real_t h2 = hx*hy;

	real_t nlx = topo.nlocal(0) - 2;
	real_t nly = topo.nlocal(1) - 2;

	real_t i1 = nlx + 1;
	real_t j1 = nly + 1;

	for (auto j : range<len_t>(1, j1)) {
		for (auto i : range<len_t>(1, i1)) {
			len_t is = igs + i;
			len_t js = jgs + j;

			real_t x = (is-1)*hx;
			real_t y = (js-1)*hy;

			b(i,j) = rhs(x, y) * h2;

		}
	}
}


static void set_solution(cedar::cdr2::mpi::grid_func & q)
{
	using namespace cedar;

	const double pi = M_PI;

	auto sol = [pi](real_t x, real_t y) {
		return sin(2*pi*x)*sin(2*pi*y);
	};

	auto & topo = q.grid();

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);

	real_t hx = 1.0 / (topo.nglobal(0) - 1);
	real_t hy = 1.0 / (topo.nglobal(1) - 1);

	for (auto j : q.range(1)) {
		for (auto i : q.range(0)) {
			len_t is = igs + i;
			len_t js = jgs + j;

			real_t x = (is-1)*hx;
			real_t y = (js-1)*hy;

			q(i,j) = sol(x,y);
		}
	}
}


TEST(MPIPoisson2, Isotropic) {
	using namespace cedar;
	using namespace cedar::cdr2;

	auto nx = 200;
	auto ny = nx;

	auto grid = util::create_topo_global(MPI_COMM_WORLD, nx, ny);

	auto so = mpi::gallery::poisson(grid);
	mpi::grid_func b(grid);

	set_problem(b,
	            [](real_t x, real_t y) {
		            const double pi = M_PI;
		            return 8*(pi*pi)*sin(2*pi*x)*sin(2*pi*y);
	            });

	auto conf = std::make_shared<config>("");
	log::init(*conf);
	mpi::solver<five_pt> bmg(so, conf);

	auto sol = bmg.solve(b);

	ASSERT_LT(std::abs(bmg.levels.template get<five_pt>(0).res.lp_norm<2>()),
	          1e-8);

	mpi::grid_func exact_sol(grid);

	set_solution(exact_sol);

	auto diff = exact_sol - sol;

	ASSERT_LT(std::abs(diff.inf_norm()),
	          1e-4);
}


TEST(MPIPoisson2, StrongX) {
	using namespace cedar;
	using namespace cedar::cdr2;

	real_t eps = 0.0001;

	auto nx = 200;
	auto ny = nx;

	auto grid = util::create_topo_global(MPI_COMM_WORLD, nx, ny);

	auto so = mpi::gallery::diag_diffusion(grid, 1, eps);
	mpi::grid_func b(grid);

	set_problem(b,
	            [eps](real_t x, real_t y) {
		            const double pi = M_PI;
		            return 4*(pi*pi)*eps*sin(2*pi*x)*sin(2*pi*y) + 4*(pi*pi)*sin(2*pi*x)*sin(2*pi*y);
	            });

	auto conf = std::make_shared<config>("");
	log::init(*conf);
	conf->set("solver.relaxation", "line-x");
	mpi::solver<five_pt> bmg(so, conf);

	auto sol = bmg.solve(b);

	ASSERT_LT(std::abs(bmg.levels.template get<five_pt>(0).res.lp_norm<2>()),
	          1e-8);

	mpi::grid_func exact_sol(grid);

	set_solution(exact_sol);

	auto diff = exact_sol - sol;

	ASSERT_LT(std::abs(diff.inf_norm()),
	          1e-4);
}


TEST(MPIPoisson2, StrongY) {
	using namespace cedar;
	using namespace cedar::cdr2;

	real_t eps = 0.0001;

	auto nx = 200;
	auto ny = nx;

	auto grid = util::create_topo_global(MPI_COMM_WORLD, nx, ny);

	auto so = mpi::gallery::diag_diffusion(grid, eps, 1);
	mpi::grid_func b(grid);

	set_problem(b,
	            [eps](real_t x, real_t y) {
		            const double pi = M_PI;
		            return 4*(pi*pi)*eps*sin(2*pi*x)*sin(2*pi*y) + 4*(pi*pi)*sin(2*pi*x)*sin(2*pi*y);
	            });

	auto conf = std::make_shared<config>("");
	log::init(*conf);
	conf->set("solver.relaxation", "line-y");
	mpi::solver<five_pt> bmg(so, conf);

	auto sol = bmg.solve(b);

	ASSERT_LT(std::abs(bmg.levels.template get<five_pt>(0).res.lp_norm<2>()),
	          1e-8);

	mpi::grid_func exact_sol(grid);

	set_solution(exact_sol);

	auto diff = exact_sol - sol;

	ASSERT_LT(std::abs(diff.inf_norm()),
	          1e-4);
}
