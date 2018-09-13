set(interface_SOURCES
  interface/c/vec.cc
  interface/c/topo.cc
  interface/c/solver.cc
  interface/c/mat.cc
  PARENT_SCOPE)

# add_library(core SHARED python/core.cc)
# set_target_properties(core PROPERTIES PREFIX "")
# target_include_directories(core PUBLIC "/usr/include/python2.7" ${CMAKE_BINARY_DIR}/include ${CMAKE_BINARY_DIR}/include/boxmg-common ${CMAKE_BINARY_DIR}/include/boxmg-2d)
# target_link_libraries(core "/usr/lib/libpython2.7.so" ${boxmg-examples_LINKER_FLAGS} boxmg-common boxmg-2d)
