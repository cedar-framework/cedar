set(sources
  relax/relax.cc
  relax/setup_relax.cc)

if(ENABLE_MPI)
  list(APPEND sources
	relax/mpi/setup_relax.cc
	relax/mpi/ml_shm.cc)
endif()

set(relax_SOURCES ${sources} PARENT_SCOPE)
