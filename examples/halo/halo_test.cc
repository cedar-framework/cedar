#include <mpi.h>
#include <memory>

#include <cedar/types.h>
#include <cedar/kernel.h>
#include <cedar/kernel_name.h>
#include <cedar/2d/kernel/mpi/factory.h>
#include <cedar/2d/kernel/mpi/registry.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/util/mpi_grid.h>

#include <cedar/util/time_log.h>


int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr2;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	timer_init(MPI_COMM_WORLD);

	config::reader conf("halo-config.json");
	log::init(conf);

	auto islocal = conf.get<bool>("grid.local", true);
	auto ndofs = conf.getvec<len_t>("grid.n");
	auto niters = conf.get<len_t>("niters");
	auto nx = ndofs[0];
	auto ny = ndofs[1];
	topo_ptr grid;
	if (islocal) {
		auto nprocs = conf.getvec<int>("grid.np");
		int npx = 0;
		int npy = 0;
		if (nprocs.size() >= 2) {
			npx = nprocs[0];
			npy = nprocs[1];
		}
		if (npx == 0 or npy == 0) {
			grid = cdr2::util::create_topo(MPI_COMM_WORLD, nx, ny);
		} else {
			int size;
			MPI_Comm_size(MPI_COMM_WORLD, &size);
			assert(size == npx*npy);
			grid = cdr2::util::create_topo(MPI_COMM_WORLD, npx, npy, nx, ny);
		}
		log::status << "Running local topo" << std::endl;
	} else {
		grid = cdr2::util::create_topo_global(MPI_COMM_WORLD, nx, ny);
		log::status << "Running global topo" << std::endl;
	}

	mpi::grid_func b(grid);
	b.set(0);

	auto kreg = std::make_shared<cedar::cdr2::kernel::mpi::registry>();

	// register alternative kernel setup/exchange kernels
	using halo_setup_t = cedar::kernel<grid_topo&, void**>;
	using halo_exchange_t = cedar::kernel<mpi::grid_func&>;
	kreg->add(kernel_name::halo_setup, "user",
	          halo_setup_t([](grid_topo & topo,
	                          void **user_ctx) -> void
	                       {
		                       timer_begin("halo-setup");
		                       // custom halo setup
		                       timer_end("halo-setup");
	                       }));

	kreg->add(kernel_name::halo_exchange, "user",
	          halo_exchange_t([](mpi::grid_func & f) -> void
	                          {
		                       timer_begin("halo-exchange");
		                       // custom halo exchange
		                       timer_end("halo-exchange");
	                       }));

	cedar::cdr2::kernel::mpi::factory::init(kreg, conf);

	void *halo_ctx;
	kreg->halo_setup(*grid, &halo_ctx);
	b.halo_ctx = halo_ctx;

	MPI_Barrier(MPI_COMM_WORLD); // synchronize before timing
	for (auto i : range<len_t>(niters)) {
		(void)i;
		kreg->halo_exchange(b);
	}

	timer_save("timings.json");

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
