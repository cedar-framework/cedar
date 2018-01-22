#include <mpi.h>
#include <iostream>

#include "ml_relax.h"

int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr2;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	timer_init(MPI_COMM_WORLD);

	config::reader conf;

	auto halo_str = conf.get<std::string>("halo-exchanger", "msg");
	if (halo_str == "msg") {
		log::status << "Using msg for halo exchange" << std::endl;
		run_test<mpi::msg_exchanger>(conf);
	} else {
		log::status << "Using tausch for halo exchange" << std::endl;
		run_test<mpi::tausch_exchanger>(conf);
	}

	timer_save("timings.json");

	// { // Output solution
	// 	int rank;
	// 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// 	std::ofstream ofile("x-" + std::to_string(rank));
	// 	ofile << x;
	// 	ofile.close();
	// }

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
