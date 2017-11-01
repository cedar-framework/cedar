#include <gtest/gtest.h>
#include <mpi.h>

#include "test_halo.h"

TEST(MPIHalo2, MSG)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	std::vector<std::array<int, 2>> procs{{ {3,3}, {3,2}, {2,3} }};

	for (auto & proc : procs) {
		int np = proc[0]*proc[1];
		int rank;
		MPI_Comm comm;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_split(MPI_COMM_WORLD, rank < np, rank, &comm);
		if (rank < np) {
			for (int per_mask = 0; per_mask < 4; per_mask++) {
				run_test<mpi::msg_exchanger>(comm, proc, per_mask);
			}
		}
	}
}


TEST(MPIHalo2, TAUSCH)
{
}
