#include <mpi.h>

#include <cedar/3d/util/topo.h>
#include <cedar/3d/mpi/gallery.h>
#include <cedar/3d/mpi/kernel_manager.h>
#include <cedar/cycle/types.h>


int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr3;

	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	config conf("openmp-relax.json");
	log::init(conf);

	auto kreg = mpi::build_kernel_manager(conf);
	auto topo = util::create_topo(conf);

	auto so = mpi::gallery::fe(topo);
	auto b = mpi::grid_func::zeros(topo);
	auto x = mpi::grid_func::ones(topo);

	// setup halo
	{
		auto & halo_service = kreg->services().get<mpi::halo_exchange>();
		std::vector<topo_ptr> topos{{so.grid_ptr()}};
		halo_service.setup(topos);
		halo_service.run(so);
	}

	relax_stencil sor(topo->nlocal(0) - 2,
	                  topo->nlocal(1) - 2,
	                  topo->nlocal(2) - 2);
	kreg->setup<mpi::point_relax>(so, sor);

	int nsweeps = conf.get<int>("nsweeps");
	for (auto i : range(nsweeps)) {
		(void)i;
		kreg->run<mpi::point_relax>(so, x, b, sor, cycle::Dir::DOWN);
	}

	MPI_Finalize();
	return 0;
}
