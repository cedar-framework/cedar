#include <gtest/gtest.h>
#include <mpi.h>

#include <boxmg/2d/mpi/gallery.h>
#include <boxmg/2d/util/topo.h>
#include <boxmg/2d/util/mpi_grid.h>
#include <boxmg/2d/mpi/solver.h>


static void set_problem(boxmg::bmg2d::mpi::grid_func & b, std::function<boxmg::real_t(boxmg::real_t,boxmg::real_t)> rhs)
{
	using namespace boxmg;
	using namespace boxmg::bmg2d;

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


TEST(MPICGSolver, Redist) {
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	auto nx = 800;
	auto ny = nx;

	auto grid0 = util::create_topo_global(MPI_COMM_WORLD, nx, ny);
	auto grid1 = util::create_topo_global(MPI_COMM_WORLD, nx, ny);

	auto so0 = mpi::gallery::poisson(grid0);
	auto so1 = mpi::gallery::poisson(grid1);
	mpi::grid_func b0(grid0);
	mpi::grid_func b1(grid1);

	auto rhsf = [](real_t x, real_t y) {
		const double pi = M_PI;
		return 8*(pi*pi)*sin(2*pi*x)*sin(2*pi*y);
	};


	set_problem(b0, rhsf);
	set_problem(b1, rhsf);

	auto conf0 = std::make_shared<config::reader>("test-cgredist-0.json");
	auto conf1 = std::make_shared<config::reader>("test-cgredist-1.json");

	log::init_level(*conf0);
	mpi::solver bmg0(std::move(so0), conf0);
	auto sol0 = bmg0.solve(b0);

	log::init_level(*conf1);
	mpi::solver bmg1(std::move(so1), conf1);
	auto sol1 = bmg1.solve(b1);

	auto diff = sol0 - sol1;

	ASSERT_LT(std::abs(diff.inf_norm()),
	          1e-10);
}
