#include <gtest/gtest.h>
#include <mpi.h>
#include <abt.h>

int main(int argc, char *argv[]) {
    int result = 0;

    ::testing::InitGoogleTest(&argc, argv);
	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
	ABT_init(argc, argv);
    result = RUN_ALL_TESTS();
	ABT_finalize();
    MPI_Finalize();

    return result;
}
