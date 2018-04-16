add_subdirectory(ftn)
add_subdirectory(util)

set(sources
  ${ftn_SOURCES}
  ${util_SOURCES}
  grid_func.cc
  stencil_op.cc
  relax_stencil.cc
  relax_planes.cc
  gallery.cc
  prolong_op.cc
  restrict_op.cc
  restrict.cc
  interp.cc
  solve_cg.cc
  kernel_manager.cc
  )

if(ENABLE_MPI)
  add_subdirectory(interface)
  list(APPEND sources
	${interface_SOURCES}
	redist/cholesky_solver.cc
	mpi/prolong_op.cc
	mpi/restrict_op.cc
	mpi/restrict.cc
	mpi/relax_planes.cc
	mpi/interp.cc
	mpi/grid_func.cc
	mpi/stencil_op.cc
    mpi/msg_exchanger.cc
	mpi/kernel_manager.cc
	mpi/setup_nog.cc
	mpi/gallery.cc)
endif()

set(3d_SOURCES ${sources} PARENT_SCOPE)
