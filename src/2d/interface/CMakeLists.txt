set(interface_SOURCES
  c/operator.cc
  c/topo.cc
  c/solver.cc
  PARENT_SCOPE)


set(interface_HEADERS
  c/operator.h
  c/topo.h
  c/solver.h
  PARENT_SCOPE)

# add_library(core SHARED python/core.cc)
# set_target_properties(core PROPERTIES PREFIX "")
# target_include_directories(core PUBLIC "/usr/include/python2.7" ${CMAKE_BINARY_DIR}/include ${CMAKE_BINARY_DIR}/include/boxmg-common ${CMAKE_BINARY_DIR}/include/boxmg-2d)
# target_link_libraries(core "/usr/lib/libpython2.7.so" ${boxmg-examples_LINKER_FLAGS} boxmg-common boxmg-2d)
