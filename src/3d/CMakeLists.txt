add_subdirectory(ftn)
add_subdirectory(inter)

set(3d_SOURCES
  ${ftn_SOURCES}
  ${inter_SOURCES}
  grid_func.cc
  grid_stencil.cc
  stencil_op.cc
  relax_stencil.cc
  residual.cc
  PARENT_SCOPE
  )
