add_unit(test-poisson-3d ${CMAKE_CURRENT_SOURCE_DIR}/test_poisson.cc ${CMAKE_CURRENT_SOURCE_DIR}/../ser_main.cc)
target_link_libraries(test-poisson-3d cedar ${cedar-deps})

if(ENABLE_MPI)
  add_subdirectory(mpi)
endif()
