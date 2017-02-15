#include <mpi.h>
#include <iostream>
#include <memory>
#include <math.h>

#include <boxmg/types.h>
#include <boxmg/2d/mpi/grid_func.h>
#include <boxmg/2d/util/topo.h>
#include <boxmg/2d/util/mpi_grid.h>
#include <boxmg/2d/mpi/solver.h>
#include <boxmg/2d/mpi/gallery.h>

#include <boxmg/util/time_log.h>


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


int main(int argc, char *argv[])
{
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	timer_init(MPI_COMM_WORLD);

	config::reader conf;

	auto islocal = conf.get<bool>("grid.local", true);
	auto ndofs = conf.getvec<len_t>("grid.n");
	auto nx = ndofs[0];
	auto ny = ndofs[1];
	topo_ptr grid;
	if (islocal) {
		auto nprocs = conf.getvec<int>("grid.np");
		int npx = 0;
		int npy = 0;
		if (nprocs.size() >= 2) {
			npx = nprocs[0];
			npy = nprocs[1];
		}
		if (npx == 0 or npy == 0) {
			grid = bmg2d::util::create_topo(MPI_COMM_WORLD, nx, ny);
		} else {
			int size;
			MPI_Comm_size(MPI_COMM_WORLD, &size);
			assert(size == npx*npy);
			grid = bmg2d::util::create_topo(MPI_COMM_WORLD, npx, npy, nx, ny);
		}
		log::status << "Running local solve" << std::endl;
	} else {
		grid = bmg2d::util::create_topo_global(MPI_COMM_WORLD, nx, ny);
		log::status << "Running global solve" << std::endl;
	}

	auto so = mpi::gallery::poisson(grid);
	mpi::grid_func b(grid);

	set_problem(b);

	mpi::solver bmg(std::move(so));

	MPI_Barrier(MPI_COMM_WORLD); // synchronize before timing solve
	auto sol = bmg.solve(b);


	mpi::grid_func exact_sol(sol.grid_ptr());
	set_solution(exact_sol);

	mpi::grid_func diff = exact_sol - sol;

	log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	timer_save("timings.json");

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
