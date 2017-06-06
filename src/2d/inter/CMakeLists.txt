set(sources
  inter/prolong_op.cc
  inter/restrict_op.cc
  inter/restrict.cc
  inter/galerkin_prod.cc
  inter/interp.cc
  inter/setup_interp.cc)

if(ENABLE_MPI)
  list(APPEND sources
	inter/mpi/prolong_op.cc
	inter/mpi/restrict_op.cc
	inter/mpi/restrict.cc
	inter/mpi/galerkin_prod.cc
	inter/mpi/interp.cc)
endif()

set(inter_SOURCES ${sources} PARENT_SCOPE)
