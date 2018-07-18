#include <cedar/2d/mpi/kernel_manager.h>
#include <cedar/2d/mpi/gallery.h>
#include <cedar/2d/util/topo.h>

int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr2;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);

	config conf("config.json");
	log::init(conf);

	auto kman = mpi::build_kernel_manager(conf);

	auto grid = util::create_topo(conf);
	auto so = mpi::gallery::poisson(grid);
	mpi::grid_func b(grid), x(grid), res(grid);
	relax_stencil sor(grid->nlocal(0) - 2, grid->nlocal(1) - 2);
	b.set(1.0);
	x.set(0.0);

	{
		auto & halo_service = kman->services().get<mpi::halo_exchange>();
		halo_service.setup({{grid}});
	}

	kman->setup<mpi::line_relax<relax_dir::x>>(so, sor);

	// for (auto i : range<std::size_t>(2))
	kman->run<mpi::line_relax<relax_dir::x>>(so, x, b, sor, res, cycle::Dir::DOWN);

	MPI_Finalize();
	return 0;
}
