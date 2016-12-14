add_par_unit(mpi-poisson-2d ${CMAKE_CURRENT_SOURCE_DIR}/test_poisson.cc ${CMAKE_CURRENT_SOURCE_DIR}/../../mpi_main.cc)
target_link_libraries(mpi-poisson-2d boxmg ${boxmg-deps})
set_par_sizes(mpi-poisson-2d 1 3 4 8)
