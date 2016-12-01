set(sources
  kernel/setup_nog.cc
  kernel/factory.cc)


if (ENABLE_MPI)
  list(APPEND sources
	kernel/mpi/factory.cc)
endif()

set(kernel_SOURCES ${sources} PARENT_SCOPE)
