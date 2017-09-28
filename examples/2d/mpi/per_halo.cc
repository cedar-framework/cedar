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

#include "per_halo.h"


static void fill_gfunc(cedar::cdr2::mpi::grid_func & b)
{
	auto & topo = b.grid();
	b.set(-1);
	for (auto j : b.range(1)) {
		for (auto i : b.range(0)) {
			b(i,j) = topo.coord(1)*topo.nproc(0) + topo.coord(0);
		}
	}
}


static void fill_stencil(cedar::cdr2::mpi::stencil_op<cedar::cdr2::five_pt> & so)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	so.set(-1);

	auto & topo = so.grid();

	for (auto j : so.range(1)) {
		for (auto i : so.range(0)) {
			so(i,j,five_pt::c) = 100*topo.coord(0) + 10*topo.coord(1);
			so(i,j,five_pt::w) = 100*topo.coord(0) + 10*topo.coord(1) + 1;
			so(i,j,five_pt::s) = 100*topo.coord(0) + 10*topo.coord(1) + 2;
		}
	}
}


int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr2;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	config::reader conf;
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
		run_test<mpi::msg_exchanger>(conf, grid, so, b);
	} else {
		log::status << "Using Tausch for halo exchange" << std::endl;
		run_test<mpi::tausch_exchanger>(conf, grid, so, b);
	}

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
