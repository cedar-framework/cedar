set(sources
  cg/solve_cg.cc)

if (ENABLE_MPI)
  list(APPEND sources
	cg/mpi/solve_cg.cc)
endif()

set(cg_SOURCES ${sources} PARENT_SCOPE)
