#include <mpi.h>
#include <iostream>
#include <memory>
#include <math.h>

#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/kernel/mpi/factory.h>
#include <cedar/2d/kernel/mpi/registry.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/util/mpi_grid.h>

#include <cedar/util/time_log.h>


static void draw(const cedar::cdr2::mpi::grid_func & b, std::ostream & os)
{
	for (auto j : b.grange(1)) {
		for (auto i : b.grange(0)) {
			if (b(i,j) < 0)
				os << '*';
			else
				os << b(i,j);
			os << " ";
		}
		os << '\n';
	}
}


int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr2;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	timer_init(MPI_COMM_WORLD);

	config::reader conf;
	log::init(conf);
	log::status << "Beginning test" << std::endl;
	auto grid = util::create_topo(conf);


	mpi::grid_func b(grid);
	auto kreg = std::make_shared<cedar::cdr2::kernel::mpi::registry>();
	cedar::cdr2::kernel::mpi::factory::init(kreg, conf);

	b.set(-1);
	void *halo_ctx;
	kreg->halo_setup(*grid, &halo_ctx);
	b.halo_ctx = halo_ctx;

	{
		for (auto j : b.range(1)) {
			for (auto i : b.range(0)) {
				b(i,j) = grid->coord(1)*grid->nproc(0) + grid->coord(0);
			}
		}

		std::ofstream ofile("output/before-" + std::to_string(grid->coord(0)) +
		                    "." + std::to_string(grid->coord(1)));
		draw(b, ofile);
		ofile.close();
	}

	kreg->halo_exchange(b);

	{
		std::ofstream ofile("output/after-" + std::to_string(grid->coord(0)) +
		                    "." + std::to_string(grid->coord(1)));
		draw(b, ofile);
		ofile.close();
	}

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
