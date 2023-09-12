add_subdirectory(ftn)

set(sources
  ${ftn_SOURCES}
  relax.cc
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
    mpi/prolong_op.cc
    mpi/restrict_op.cc
    mpi/gallery.cc
    mpi/msg_exchanger.cc
    mpi/setup_nog.cc
    mpi/restrict.cc
    mpi/interp.cc
    mpi/kernel_manager.cc
    mpi/redist_solver.cc
    mpi/solver.cc
    redist/cholesky_solver.cc
    mpi/tausch_exchanger.cc)
endif()

list(APPEND sources
  gpu/grid_func.cc
  gpu/stencil_op.cc
  gpu/prolong_op.cc
  gpu/restrict_op.cc
  gpu/setup_nog.cc
  gpu/restrict.cc
  gpu/interp.cc
  gpu/kernel_manager.cc
  gpu/redist_solver.cc
  gpu/solver.cc
  gpu/gallery.cc
  gpu/tausch_exchanger.cc
  gpu/msg_exchanger.cc)

set(2d_SOURCES ${sources} PARENT_SCOPE)
set(2d_FTL_SOURCES ${ftn_FTL_SOURCES} PARENT_SCOPE)
