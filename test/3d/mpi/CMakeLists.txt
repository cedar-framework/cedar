add_par_unit(test-mpi-poisson-3d ${CMAKE_CURRENT_SOURCE_DIR}/test_poisson.cc ${CMAKE_CURRENT_SOURCE_DIR}/../../mpi_main.cc)
target_link_libraries(test-mpi-poisson-3d cedar ${cedar-deps})
set_par_sizes(test-mpi-poisson-3d 1 3 4 8 9)


add_par_unit(test-mpi-relax-3d ${CMAKE_CURRENT_SOURCE_DIR}/test_relax.cc ${CMAKE_CURRENT_SOURCE_DIR}/../../mpi_main.cc)
target_link_libraries(test-mpi-relax-3d cedar ${cedar-deps})
set_par_sizes(test-mpi-relax-3d 1 3 8 9)

add_par_unit(test-mpi-cgredist-3d ${CMAKE_CURRENT_SOURCE_DIR}/test_cgredist.cc ${CMAKE_CURRENT_SOURCE_DIR}/../../mpi_main.cc)
target_link_libraries(test-mpi-cgredist-3d cedar ${cedar-deps})
set_par_sizes(test-mpi-cgredist-3d 8)

configure_file(test-cgredist-0.json ${CMAKE_BINARY_DIR}/test-cgredist3-0.json COPYONLY)
configure_file(test-cgredist-1.json ${CMAKE_BINARY_DIR}/test-cgredist3-1.json COPYONLY)
