add_subdirectory(ftn)
add_subdirectory(inter)
add_subdirectory(relax)
add_subdirectory(kernel)
add_subdirectory(cg)
add_subdirectory(util)
add_subdirectory(interface)

set(sources
  ${ftn_SOURCES}
  ${inter_SOURCES}
  ${relax_SOURCES}
  ${kernel_SOURCES}
  ${cg_SOURCES}
  ${util_SOURCES}
  ${interface_SOURCES}
  grid_func.cc
  grid_stencil.cc
  stencil_op.cc
  relax_stencil.cc
  residual.cc
  solver.cc
  gallery.cc
  )

if(ENABLE_MPI)
  list(APPEND sources
	mpi/grid_func.cc
	mpi/stencil_op.cc
	mpi/halo.cc
	mpi/solver.cc
	mpi/residual.cc
	mpi/gallery.cc
	mpi/redist_solver.cc)
endif()

set(3d_SOURCES ${sources} PARENT_SCOPE)
