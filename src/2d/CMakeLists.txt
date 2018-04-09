add_subdirectory(ftn)

set(sources
  ${ftn_SOURCES}
  restrict_op.cc
  prolong_op.cc
  restrict.cc
  interp.cc
  grid_func.cc
  relax_stencil.cc
  stencil_op.cc
  gallery.cc
  solve_cg.cc
  kernel_manager.cc
  )

if(ENABLE_MPI)
  add_subdirectory(interface)
  add_subdirectory(util)

  list(APPEND sources
	${interface_SOURCES}
	${util_SOURCES}
	mpi/grid_func.cc
	mpi/stencil_op.cc
	mpi/gallery.cc
	mpi/msg_exchanger.cc
	mpi/tausch_exchanger.cc
	redist/cholesky_solver.cc
	)
endif()

set(2d_SOURCES ${sources} PARENT_SCOPE)
