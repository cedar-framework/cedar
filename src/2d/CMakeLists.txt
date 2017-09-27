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
  relax_stencil.cc
  residual.cc
  stencil_op.cc
  gallery.cc
  )

if(ENABLE_MPI)
  add_subdirectory(interface)
  add_subdirectory(util)

  list(APPEND sources
	${interface_SOURCES}
	${util_SOURCES}
	mpi/grid_func.cc
	mpi/stencil_op.cc
	mpi/redist_solver.cc
	mpi/gallery.cc
	mpi/msg_exchanger.cc
	mpi/tausch_exchanger.cc
	)
endif()

set(2d_SOURCES ${sources} PARENT_SCOPE)
