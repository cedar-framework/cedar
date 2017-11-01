#include <gtest/gtest.h>
#include <mpi.h>

#include "test_halo.h"

TEST(MPIHalo2, MSG)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	std::vector<std::array<int, 2>> procs{{ {3,3}, {3,2}, {2,3} }};
	test_driver<mpi::msg_exchanger>(procs);
}


TEST(MPIHalo2, TAUSCH)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	std::vector<std::array<int, 2>> procs{{ {3,3}, {3,2}, {2,3} }};
	test_driver<mpi::tausch_exchanger>(procs);
}
