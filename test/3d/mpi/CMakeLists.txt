add_par_unit(mpi-poisson-3d ${CMAKE_CURRENT_SOURCE_DIR}/test_poisson.cc ${CMAKE_CURRENT_SOURCE_DIR}/../../mpi_main.cc)
target_link_libraries(mpi-poisson-3d boxmg ${boxmg-deps})
set_par_sizes(mpi-poisson-3d 1 3 4 8 9)


add_par_unit(mpi-relax-3d ${CMAKE_CURRENT_SOURCE_DIR}/test_relax.cc ${CMAKE_CURRENT_SOURCE_DIR}/../../mpi_main.cc)
target_link_libraries(mpi-relax-3d boxmg ${boxmg-deps})
set_par_sizes(mpi-relax-3d 1)

