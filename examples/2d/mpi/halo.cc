#include <mpi.h>
#include <iostream>
#include <memory>
#include <math.h>

#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/kernel/mpi/registry.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/mpi/tausch_exchanger.h>
#include <cedar/2d/mpi/msg_exchanger.h>
#include <cedar/kernel_params.h>

#include <cedar/util/time_log.h>

#include "halo.h"


int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr2;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	timer_init(MPI_COMM_WORLD);

	config::reader conf("halo-config.json");
	log::init(conf);
	log::status << "Beginning test" << std::endl;
	auto grid = util::create_topo(conf);

	mpi::grid_func b(grid);
	mpi::stencil_op<five_pt> so(grid);

	fill_gfunc(b);
	fill_stencil(so);

	auto exchanger_str = conf.get<std::string>("halo-exchanger", "msg");

	if (exchanger_str == "msg") {
		log::status << "Using MSG for halo exchange" << std::endl;
		run_halo<mpi::msg_exchanger>(conf, grid, so, b);
	} else {
		log::status << "Using Tausch for halo exchange" << std::endl;
		run_halo<mpi::tausch_exchanger>(conf, grid, so, b);
	}

	timer_save("timings-" + exchanger_str + ".json");

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
