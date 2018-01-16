#include <mpi.h>
#include <iostream>
#include <memory>
#include <math.h>

#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/mpi/gallery.h>

#include <cedar/util/time_log.h>


static void set_problem(cedar::cdr2::mpi::grid_func & b)
{
	using namespace cedar;
	using namespace cedar::cdr2;

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


int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr2;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	timer_init(MPI_COMM_WORLD);

	config::reader conf;

	auto grid = util::create_topo(conf);
	grid->grow(2);

	auto niters = conf.get<std::size_t>("niters");

	auto so = mpi::gallery::poisson(grid);
	mpi::grid_func b(grid), x(grid), res(grid);
	relax_stencil sor(grid->nlocal(0) - 2, grid->nlocal(1) - 2);

	set_problem(b);

	kernel::mpi::registry<mpi::msg_exchanger> kreg(conf);

	{ // setup halo
		std::vector<topo_ptr> topos;
		topos.push_back(so.grid_ptr());
		kreg.halo_setup(topos);
		kreg.halo_stencil_exchange(so);
	}

	kreg.setup_relax_x(so, sor);

	timer_begin("line-relax");
	for (auto i : range(niters)) {
		kreg.relax_lines_x(so, x, b, sor, res, cycle::Dir::DOWN);
	}
	timer_end("line-relax");

	MPI_Barrier(MPI_COMM_WORLD); // synchronize before timing solve

	timer_save("timings.json");


	// { // Output solution
	// 	int rank;
	// 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// 	std::ofstream ofile("x-" + std::to_string(rank));
	// 	ofile << x;
	// 	ofile.close();
	// }

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
