#include <mpi.h>
#include <iostream>
#include <memory>
#include <chrono>
#include <thread>
#include <math.h>

#include <boxmg/types.h>
#include <boxmg/kernel.h>
#include <boxmg/kernel_name.h>
#include <boxmg/2d/mpi/grid_func.h>
#include <boxmg/2d/util/topo.h>
#include <boxmg/2d/util/mpi_grid.h>
#include <boxmg/2d/mpi/solver.h>

#include <boxmg/util/time_log.h>

int main(int argc, char *argv[])
{
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	timer_init(MPI_COMM_WORLD);

	config::reader conf;

	auto grid = bmg2d::util::create_topo_global(MPI_COMM_WORLD, 100, 100);

	auto so = mpi::stencil_op(grid);

	mpi::grid_func b(so.grid_ptr());

	int rank;
	MPI_Comm_rank(so.grid().comm, &rank);

	mpi::solver bmg(std::move(so));

	auto kreg = bmg.kernel_registry();

	using rkern_t = boxmg::kernel<const mpi::stencil_op&,
	                              mpi::grid_func &,
	                              const mpi::grid_func&,
	                              const bmg2d::relax_stencil&,
	                              cycle::Dir>;

	kreg->add(kernel_name::relax, "user", rkern_t([](const mpi::stencil_op&,
	                                                 mpi::grid_func&,
	                                                 const mpi::grid_func&,
	                                                 const bmg2d::relax_stencil&,
	                                                 cycle::Dir) -> void {
		                                              timer_begin("user");
		                                              std::this_thread::sleep_for(std::chrono::seconds(1));
		                                              timer_end("user");
	                                              }));
	kreg->set(kernel_name::relax, "user");

	auto sol = bmg.solve(b);

	timer_save("timings.json");

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();

	return 0;
}
