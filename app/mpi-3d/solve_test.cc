#include <mpi.h>
#include <iostream>

#include <boxmg/types.h>
#include <boxmg/3d/mpi/grid_func.h>
#include <boxmg/3d/util/topo.h>
#include <boxmg/3d/mpi/solver.h>
#include <boxmg/3d/mpi/gallery.h>


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


static void set_solution(boxmg::bmg3::mpi::grid_func & q)
{
	using namespace boxmg;

	const double pi = M_PI;

	auto sol = [pi](real_t x, real_t y, real_t z) {
		return sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z);
	};

	auto & topo = q.grid();

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);
	real_t kgs = topo.is(2);

	real_t hx = 1.0 / (topo.nglobal(0) - 1);
	real_t hy = 1.0 / (topo.nglobal(1) - 1);
	real_t hz = 1.0 / (topo.nglobal(2) - 1);

	for (auto k : q.range(2)) {
		for (auto j : q.range(1)) {
			for (auto i : q.range(0)) {
				len_t is = igs + i;
				len_t js = jgs + j;
				len_t ks = kgs + k;

				real_t x = (is-1)*hx;
				real_t y = (js-1)*hy;
				real_t z = (ks-1)*hz;

				q(i,j,k) = sol(x,y,z);
			}
		}
	}
}


int main(int argc, char *argv[])
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

	int provided, rank;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	timer_init(MPI_COMM_WORLD);
	config::reader conf("config.json");

	auto islocal = conf.get<bool>("grid.local", true);
	auto ndofs = conf.getvec<len_t>("grid.n");
	auto nx = ndofs[0];
	auto ny = ndofs[1];
	auto nz = ndofs[2];
	topo_ptr grid;
	if (islocal) {
		auto np = conf.getvec<int>("grid.np");
		if (np.size() >= 3) {
			int size;
			MPI_Comm_size(MPI_COMM_WORLD, &size);
			assert(size == np[0]*np[1]*np[2]);
			grid = bmg3::util::create_topo(MPI_COMM_WORLD, np[0], np[1], np[2],
			                              nx, ny, nz);
		} else {
			grid = bmg3::util::create_topo(MPI_COMM_WORLD, nx, ny, nz);
		}

		log::status << "Running local solve" << std::endl;
	} else {
		grid = bmg3::util::create_topo_global(MPI_COMM_WORLD, nx, ny, nz);
		log::status << "Running global solve" << std::endl;
	}

	auto so = mpi::gallery::poisson(grid);
	mpi::grid_func b(grid);

	set_problem(b);
	mpi::solver bmg(std::move(so));

	MPI_Barrier(MPI_COMM_WORLD); // synchronize before timing solve
	auto x = bmg.solve(b);

	mpi::grid_func exact_sol(x.grid_ptr());

	set_solution(exact_sol);

	mpi::grid_func diff = exact_sol - x;

	log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	// {
	// 	std::ofstream ofile;
	// 	auto & topo = bmg.level(-1).A.grid();
	// 	ofile.open("op-" + std::to_string(topo.coord(0)+1) + "." +
	// 	           std::to_string(topo.coord(1)+1) + "." + std::to_string(topo.coord(2)+1) + ".txt", std::ios::out | std::ios::trunc | std::ios::binary);
	// 	ofile << bmg.level(0).A;
	// 	ofile.close();
	// }

	timer_save("timings.json");
	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
