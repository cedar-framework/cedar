#include <gtest/gtest.h>
#include <mpi.h>

#include <boxmg/2d/mpi/gallery.h>
#include <boxmg/2d/util/topo.h>
#include <boxmg/2d/util/mpi_grid.h>
#include <boxmg/2d/mpi/solver.h>


static void set_problem(boxmg::bmg2d::mpi::grid_func & b)
{
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	const double pi = M_PI;

	auto rhs = [pi](real_t x, real_t y) {
		return 8*(pi*pi)*sin(2*pi*x)*sin(2*pi*y);
	};

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


static void set_solution(boxmg::bmg2d::mpi::grid_func & q)
{
	using namespace boxmg;

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
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	auto nx = 200;
	auto ny = nx;

	auto grid = util::create_topo_global(MPI_COMM_WORLD, nx, ny);

	auto so = mpi::gallery::poisson(grid);
	mpi::grid_func b(grid);

	set_problem(b);

	auto conf = std::make_shared<config::reader>("");
	mpi::solver bmg(std::move(so), conf);

	auto sol = bmg.solve(b);

	ASSERT_LT(std::abs(bmg.level(-1).res.lp_norm<2>()),
	          1e-8);

	mpi::grid_func exact_sol(grid);

	set_solution(exact_sol);

	auto diff = exact_sol - sol;

	ASSERT_LT(std::abs(diff.inf_norm()),
	          1e-4);
}
