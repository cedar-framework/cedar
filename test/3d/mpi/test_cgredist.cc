#include <gtest/gtest.h>
#include <mpi.h>

#include <boxmg/types.h>
#include <boxmg/3d/mpi/grid_func.h>
#include <boxmg/3d/util/topo.h>
#include <boxmg/3d/mpi/gallery.h>
#include <boxmg/3d/mpi/solver.h>


static void set_problem(boxmg::bmg3::mpi::grid_func & b)
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

	const double pi = M_PI;

	auto rhs = [pi](real_t x, real_t y, real_t z) {
		return 12*(pi*pi)*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z);
	};

	auto & topo = b.grid();

	b.set(0);

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);
	real_t kgs = topo.is(2);

	real_t hx = 1.0 / (topo.nglobal(0) - 1);
	real_t hy = 1.0 / (topo.nglobal(1) - 1);
	real_t hz = 1.0 / (topo.nglobal(2) - 1);

	real_t h2 = hx*hy*hz;

	real_t nlx = topo.nlocal(0) - 2;
	real_t nly = topo.nlocal(1) - 2;
	real_t nlz = topo.nlocal(2) - 2;

	real_t i1 = nlx + 1;
	real_t j1 = nly + 1;
	real_t k1 = nlz + 1;

	for (auto k : range<len_t>(1, k1)) {
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				len_t is = igs + i;
				len_t js = jgs + j;
				len_t ks = kgs + k;

				real_t x = (is-1)*hx;
				real_t y = (js-1)*hy;
				real_t z = (ks-1)*hz;

				b(i,j,k) = rhs(x, y, z) * h2;
			}
		}
	}
}


TEST(MPICGSolver3, Redist) {
	using namespace boxmg;
	using namespace boxmg::bmg3;

	auto nx = 200;
	auto ny = nx;
	auto nz = nx;

	auto grid0 = util::create_topo_global(MPI_COMM_WORLD, nx, ny, nz);
	auto grid1 = util::create_topo_global(MPI_COMM_WORLD, nx, ny, nz);

	auto so0 = mpi::gallery::poisson(grid0);
	auto so1 = mpi::gallery::poisson(grid1);
	mpi::grid_func b0(grid0);
	mpi::grid_func b1(grid1);

	set_problem(b0);
	set_problem(b1);

	auto conf0 = std::make_shared<config::reader>("test-cgredist3-0.json");
	auto conf1 = std::make_shared<config::reader>("test-cgredist3-1.json");

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
