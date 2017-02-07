set(sources
  inter/prolong_op.cc
  inter/restrict_op.cc
  inter/setup_interp.cc
  inter/galerkin_prod.cc
  inter/restrict.cc
  inter/interp.cc)

if(ENABLE_MPI)
  list(APPEND sources
	inter/mpi/setup_interp.cc
	inter/mpi/galerkin_prod.cc
	inter/mpi/restrict.cc
	inter/mpi/interp.cc
	inter/mpi/prolong_op.cc
	inter/mpi/restrict_op.cc)
endif()

set(inter_SOURCES ${sources} PARENT_SCOPE)
