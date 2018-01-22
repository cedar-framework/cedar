#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/mpi/gallery.h>
#include <cedar/util/time_log.h>
#include <cedar/2d/mpi/tausch_exchanger.h>
#include <cedar/2d/util/topo.h>
#include <cedar/util/time_log.h>


namespace cedar { namespace cdr2 {

template<class halo>
void run_test(config::reader & conf)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	auto grid = util::create_topo(conf);

	auto niters = conf.get<std::size_t>("niters");

	auto so = mpi::gallery::poisson(grid);
	mpi::grid_func b(grid), x(grid), res(grid);
	relax_stencil sor(grid->nlocal(0) - 2, grid->nlocal(1) - 2);

	b.set(0);
	x.set(1);

	kernel::mpi::registry<halo> kreg(conf);

	{ // setup halo
		std::vector<topo_ptr> topos;
		topos.push_back(so.grid_ptr());
		kreg.halo_setup(topos);
		kreg.halo_stencil_exchange(so);
	}

	kreg.setup_relax_x(so, sor);

	MPI_Barrier(MPI_COMM_WORLD); // synchronize before timing solve

	timer_begin("line-relax");
	for (auto i : range(niters)) {
		kreg.relax_lines_x(so, x, b, sor, res, cycle::Dir::DOWN);
	}
	timer_end("line-relax");
}

}}
