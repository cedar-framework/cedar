#include <cedar/3d/mpi/kernel_manager.h>
#include <cedar/3d/mpi/gallery.h>
#include <cedar/3d/util/topo.h>

int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr3;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);

	config conf("config.json");
	log::init(conf);

	auto kman = mpi::build_kernel_manager(conf);

	auto grid = util::create_topo(conf);
	auto so = mpi::gallery::poisson(grid);
	mpi::grid_func b(grid), x(grid);
	b.set(1.0);
	x.set(0.0);

	{
		auto & halo_service = kman->services().get<mpi::halo_exchange>();
		halo_service.setup({{grid}});
	}

	kman->setup<mpi::plane_relax<relax_dir::xy>>(so);

	// for (auto i : range<std::size_t>(2))
	kman->run<mpi::plane_relax<relax_dir::xy>>(so, x, b, cycle::Dir::DOWN);

	MPI_Finalize();
	return 0;
}
