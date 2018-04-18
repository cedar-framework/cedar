add_unit(test-poisson-3d ${CMAKE_CURRENT_SOURCE_DIR}/test_poisson.cc ${CMAKE_CURRENT_SOURCE_DIR}/../ser_main.cc)
target_link_libraries(test-poisson-3d cedar)

add_unit(test-planes-3d ${CMAKE_CURRENT_SOURCE_DIR}/pyplanes.c ${CMAKE_CURRENT_SOURCE_DIR}/test_planes.cc)
target_link_libraries(test-planes-3d cedar)

configure_file(test-planes-ser.json ${CMAKE_BINARY_DIR}/test-planes-ser.json COPYONLY)

if(ENABLE_MPI)
  add_subdirectory(mpi)
endif()
