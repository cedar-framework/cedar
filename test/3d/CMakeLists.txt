add_unit(poisson-3d ${CMAKE_CURRENT_SOURCE_DIR}/test_poisson.cc ${CMAKE_CURRENT_SOURCE_DIR}/../ser_main.cc)
target_link_libraries(poisson-3d boxmg ${boxmg-deps})

if(ENABLE_MPI)
  add_subdirectory(mpi)
endif()
