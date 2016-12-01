add_subdirectory(cg)
add_subdirectory(ftn)
add_subdirectory(inter)
add_subdirectory(kernel)
add_subdirectory(relax)



set(sources
  ${cg_SOURCES}
  ${ftn_SOURCES}
  ${inter_SOURCES}
  ${kernel_SOURCES}
  ${relax_SOURCES}
  grid_func.cc
  grid_stencil.cc
  stencil_op.cc
  boundary_iterator.cc
  matvec.cc
  relax_stencil.cc
  residual.cc
  solver.cc
  stencil_op.cc
  )

if(ENABLE_MPI)
  add_subdirectory(interface)
  add_subdirectory(util)

  list(APPEND sources
	${interface_SOURCES}
	${util_SOURCES}
	mpi/grid_func.cc
	mpi/halo.cc
	mpi/solver.cc
	mpi/stencil_op.cc
	mpi/redist_solver.cc
	mpi/residual.cc
	)
endif()

set(2d_SOURCES ${sources} PARENT_SCOPE)
