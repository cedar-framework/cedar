#include <mpi.h>
#include <abt.h>

#include <iostream>

#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/mpi/gallery.h>
#include <cedar/2d/mpi/kernel_manager.h>

#include <cedar/util/time_log.h>


#include "coll.h"

using namespace cedar;
using namespace cedar::cdr2;

int main(int argc, char *argv[])
{
	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);

	timer_init(MPI_COMM_WORLD);

	auto conf = std::make_shared<config>("config.json");
	log::init(*conf);

	auto nplanes = conf->get<int>("nplanes");

	auto grid = util::create_topo(*conf);
	std::vector<mpi::kman_ptr> kmans;
	kmans.reserve(nplanes);
	std::vector<int> split_keys;

	auto registration = [&](service_manager<mpi::stypes> & services) {
		services.add<mpi::message_passing, master_splitter>("master", &split_keys);
		services.add<mpi::message_passing, worker_splitter>("worker", &split_keys);
	};

	for (auto i : range(nplanes)) {
		(void)i;
		kmans.push_back(mpi::build_kernel_manager(*conf));

		auto kman = kmans.back();
		auto & sman = kman->services();
		sman.set_user_reg(registration);
	}

	kmans[0]->services().set<mpi::message_passing>("master");
	for (auto i : range(1, nplanes))
		kmans[i]->services().set<mpi::message_passing>("worker");

	MPI_Comm *newcomm;
	newcomm = new MPI_Comm[2*nplanes];
	for (auto i : range(nplanes)) {
		auto & mp = kmans[i]->services().get<mpi::message_passing>();
		mp.comm_split(grid->comm, grid->coord(0), grid->coord(1), &newcomm[i*2]);
		mp.comm_split(grid->comm, grid->coord(1), grid->coord(0), &newcomm[i*2+1]);

		if (i == 0)
			printf("master: %p %p\n", newcomm[i*2], newcomm[i*2+1]);
		else if (i == 1)
			printf("worker: %p %p\n", newcomm[i*2], newcomm[i*2+1]);
	}

	{
		auto newkman = mpi::build_kernel_manager(*conf);
		auto & oldsman = kmans[1]->services();
		auto & newsman = newkman->services();
		newsman.set_user_reg(oldsman.get_user_reg());
		newsman.set<mpi::message_passing>(oldsman.get_key<mpi::message_passing>());

		MPI_Comm comm;
		auto & mp = newsman.get<mpi::message_passing>();
		mp.comm_split(grid->comm, grid->coord(0), grid->coord(1), &comm);
		printf("worker: %p\n", comm);
	}

	delete[] newcomm;

	MPI_Finalize();

	return 0;
}
