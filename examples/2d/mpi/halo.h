#ifndef CEDAR_EXAMPLES_2D_MPI_HALO_H
#define CEDAR_EXAMPLES_2D_MPI_HALO_H

#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/kernel/mpi/registry.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/mpi/tausch_exchanger.h>
#include <cedar/2d/mpi/msg_exchanger.h>
#include <cedar/kernel_params.h>

#include <cedar/util/time_log.h>

namespace cedar
{
	using namespace cedar::cdr2;

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


template<class sten>
static void fill_stencil(cedar::cdr2::mpi::stencil_op<sten> & so)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	so.set(-1);

	auto & topo = so.grid();

	for (auto j : so.range(1)) {
		for (auto i : so.range(0)) {
			for (int k = 0; k < stencil_ndirs<sten>::value; k++) {
				so(i,j,static_cast<sten>(k)) = 100*topo.coord(0) + 10*topo.coord(1) + k;
			}
		}
	}
}


template<class halo_exchanger>
	void run_halo(config::reader & conf,
	              std::shared_ptr<grid_topo> grid,
	              mpi::stencil_op<five_pt> & so,
	              mpi::grid_func & b)
{
	auto params = build_kernel_params(conf);

	halo_exchanger halo(*params, {grid});

	auto niters = conf.get<int>("niters", 100);
	log::status << "# iterations: " << niters << std::endl;

	timer_begin("halo-func");
	for (auto i : range<int>(niters)) {
		(void)i;
		halo.exchange(b);
	}
	timer_end("halo-func");

	timer_begin("halo-sten");
	for (auto i : range<int>(niters)) {
		(void)i;
		halo.exchange(so);
	}
	timer_end("halo-sten");
}

}
#endif
