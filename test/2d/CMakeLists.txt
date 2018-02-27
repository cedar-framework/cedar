add_unit(test-relax-2d ${CMAKE_CURRENT_SOURCE_DIR}/pyrelax.c ${CMAKE_CURRENT_SOURCE_DIR}/test_relax.cc)
target_link_libraries(test-relax-2d cedar ${cedar-deps})
add_unit(test-poisson-2d ${CMAKE_CURRENT_SOURCE_DIR}/test_poisson.cc ${CMAKE_CURRENT_SOURCE_DIR}/../ser_main.cc)
target_link_libraries(test-poisson-2d cedar ${cedar-deps})

if(ENABLE_MPI)
  add_subdirectory(mpi)
endif()
