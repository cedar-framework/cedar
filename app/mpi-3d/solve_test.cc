#include <mpi.h>
#include <iostream>

#include <boxmg/types.h>
#include <boxmg/3d/mpi/grid_func.h>
#include <boxmg/3d/util/topo.h>
#include <boxmg/3d/mpi/solver.h>

extern "C" {
	using namespace boxmg;
	void putf(real_t *so, real_t *qf,
	          len_t nlx, len_t nly, len_t lnz,
	          len_t ngx, len_t ngy, len_t ngz,
	          len_t igs, len_t jgs, len_t kgs,
	          real_t hx, real_t hy, real_t hz);
	void put_sol(real_t *q,
	             len_t nlx, len_t nly, len_t nlz,
	             len_t igs, len_t jgs, len_t kgs,
	             real_t hx, real_t hy, real_t hz);
}

static void set_problem(boxmg::bmg3::mpi::stencil_op & so,
                        boxmg::bmg3::mpi::grid_func & b)
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

	real_t hx = (1.0/(so.grid().nglobal(0)-1));
	real_t hy = (1.0/(so.grid().nglobal(1)-1));
	real_t hz = (1.0/(so.grid().nglobal(2)-1));

	auto & sten = so.stencil();
	sten.five_pt() = true;
	auto & topo = so.grid();

	putf(so.data(), b.data(),
	     topo.nlocal(0)-2, topo.nlocal(1)-2, topo.nlocal(2)-2,
	     topo.nglobal(0)-2, topo.nglobal(1)-2, topo.nglobal(2)-2,
	     topo.is(0), topo.is(1), topo.is(2),
	     hx, hy, hz);
}


static void set_solution(boxmg::bmg3::mpi::grid_func & x)
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

	auto & topo = x.grid();

	real_t hx = (1.0/(topo.nglobal(0)-1));
	real_t hy = (1.0/(topo.nglobal(1)-1));
	real_t hz = (1.0/(topo.nglobal(2)-1));

	put_sol(x.data(),
	        topo.nlocal(0)-2, topo.nlocal(1)-2, topo.nlocal(2)-2,
	        topo.is(0), topo.is(1), topo.is(2),
	        hx, hy, hz);
}

int main(int argc, char *argv[])
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

	int provided, rank;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	timer_init(MPI_COMM_WORLD);
	config::reader conf;

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

	auto so = mpi::stencil_op(grid);
	mpi::grid_func b(so.grid_ptr());

	set_problem(so, b);
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
