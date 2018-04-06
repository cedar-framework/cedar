add_subdirectory(ftn)
add_subdirectory(inter)
#add_subdirectory(relax)
add_subdirectory(kernel)
add_subdirectory(cg)
add_subdirectory(util)

set(sources
  ${ftn_SOURCES}
  ${inter_SOURCES}
  #${relax_SOURCES}
  ${kernel_SOURCES}
  ${cg_SOURCES}
  ${util_SOURCES}
  grid_func.cc
  stencil_op.cc
  relax_stencil.cc
  solver.cc
  gallery.cc
  )

if(ENABLE_MPI)
  add_subdirectory(interface)
  list(APPEND sources
	${interface_SOURCES}
	redist/cholesky_solver.cc
	mpi/grid_func.cc
	mpi/stencil_op.cc
    mpi/msg_exchanger.cc
	mpi/gallery.cc)
endif()

set(3d_SOURCES ${sources} PARENT_SCOPE)
