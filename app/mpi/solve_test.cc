#include <mpi.h>
#include <iostream>
#include <memory>
#include <math.h>

#include <boxmg/types.h>
#include <boxmg/2d/mpi/grid_func.h>
#include <boxmg/2d/util/topo.h>
#include <boxmg/2d/util/mpi_grid.h>
#include <boxmg/2d/mpi/solver.h>

#include <boxmg/util/time_log.h>

extern "C" {
	using namespace boxmg;
	void putf(real_t *so, real_t *qf,
	          len_t nlx, len_t nly,
	          len_t ngx, len_t ngy,
	          len_t igs, len_t jgs,
	          real_t hx, real_t hy);
}


static void set_problem(boxmg::bmg2d::mpi::stencil_op & so,
                        boxmg::bmg2d::mpi::grid_func & b)
{
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	real_t hx = (1.0/(so.grid().nglobal(0)-1));
	real_t hy = (1.0/(so.grid().nglobal(1)-1));

	auto & sten = so.stencil();
	sten.five_pt() = true;
	auto & topo = so.grid();

	b.set(0);
	sten.set(0);

	putf(so.data(), b.data(),
	     topo.nlocal(0) - 2, topo.nlocal(1) - 2,
	     topo.nglobal(0) - 2, topo.nglobal(1) - 2,
	     topo.is(0), topo.is(1),
	     hx, hy);
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

	auto so = mpi::stencil_op(grid);

	const double pi = M_PI;

	real_t hx = (1.0/(so.grid().nglobal(0)-1));
	real_t hy = (1.0/(so.grid().nglobal(1)-1));
	real_t h2 = hx * hy;

	mpi::grid_func b(so.grid_ptr());

	set_problem(so, b);

	int rank;
	MPI_Comm_rank(so.grid().comm, &rank);

	mpi::solver bmg(std::move(so));

	auto sol = bmg.solve(b);


	mpi::grid_func exact_sol(sol.grid_ptr());

	real_t y = sol.grid().is(1)*hy;
	for (auto j: exact_sol.range(1)) {
		real_t x = sol.grid().is(0)*hx;
		for (auto i: exact_sol.range(0)) {
			exact_sol(i,j) = sin(2*pi*x) * sin(2*pi*y);
			x += hx;
		}
		y+= hy;
	}

	mpi::grid_func diff = exact_sol - sol;

	log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	timer_save("timings.json");

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
