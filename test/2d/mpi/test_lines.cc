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

TEST(MPILines, MLRelaxX)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	real_t eps = 0.0001;

	auto nx = 800;
	auto ny = 10;

	auto grid = util::create_topo_global(MPI_COMM_WORLD, nx, ny);

	auto so = mpi::gallery::diag_diffusion(grid, 1, eps);


	mpi::grid_func b(grid);

	set_problem(b,
	            [eps](real_t x, real_t y) {
		            const double pi = M_PI;
		            return 4*(pi*pi)*sin(2*pi*x)*sin(2*pi*y) + 4*(pi*pi)*eps*sin(2*pi*x)*sin(2*pi*y);
	            });

	auto conf2 = std::make_shared<config>("");
	log::init(*conf2);
	conf2->set("solver.relaxation", "line-x");
	mpi::solver<five_pt> bmg2(so, conf2);

	auto confn = std::make_shared<config>("");
	confn->set("solver.relaxation", "line-x");
	confn->set("solver.ml-relax.enabled", true);
	mpi::solver<five_pt> bmgn(so, confn);

	auto sol2 = bmg2.solve(b);
	auto soln = bmgn.solve(b);

	auto diff = sol2 - soln;

	ASSERT_LT(std::abs(diff.inf_norm()), 1e-10);
}


TEST(MPILines, MLRelaxY)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	real_t eps = 0.0001;

	auto nx = 10;
	auto ny = 800;

	auto grid = util::create_topo_global(MPI_COMM_WORLD, nx, ny);

	auto so = mpi::gallery::diag_diffusion(grid, eps, 1);


	mpi::grid_func b(grid);

	set_problem(b,
	            [eps](real_t x, real_t y) {
		            const double pi = M_PI;
		            return 4*(pi*pi)*eps*sin(2*pi*x)*sin(2*pi*y) + 4*(pi*pi)*sin(2*pi*x)*sin(2*pi*y);
	            });

	auto conf2 = std::make_shared<config>("");
	log::init(*conf2);
	conf2->set("solver.relaxation", "line-y");
	mpi::solver<five_pt> bmg2(so, conf2);

	auto confn = std::make_shared<config>("");
	confn->set("solver.relaxation", "line-y");
	confn->set("solver.ml-relax.enabled", true);
	mpi::solver<five_pt> bmgn(so, confn);

	auto sol2 = bmg2.solve(b);
	auto soln = bmgn.solve(b);

	auto diff = sol2 - soln;

	ASSERT_LT(std::abs(diff.inf_norm()), 1e-10);
}
