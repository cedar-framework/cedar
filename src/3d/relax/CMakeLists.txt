set(sources
  relax/setup_relax.cc
  relax/setup_planes.cc
  relax/relax.cc
  relax/relax_planes.cc)

if(ENABLE_MPI)
  list(APPEND sources
	relax/mpi/relax.cc
	relax/mpi/setup_relax.cc)
endif()

set(relax_SOURCES ${sources} PARENT_SCOPE)
