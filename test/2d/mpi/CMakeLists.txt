add_par_unit(test-mpi-poisson-2d ${CMAKE_CURRENT_SOURCE_DIR}/test_poisson.cc ${CMAKE_CURRENT_SOURCE_DIR}/../../mpi_main.cc)
target_link_libraries(test-mpi-poisson-2d cedar ${cedar-deps})
set_par_sizes(test-mpi-poisson-2d 1 3 4 8)

add_par_unit(test-mpi-cgsolve-2d ${CMAKE_CURRENT_SOURCE_DIR}/test_cgsolve.cc ${CMAKE_CURRENT_SOURCE_DIR}/../../mpi_main.cc)
target_link_libraries(test-mpi-cgsolve-2d cedar ${cedar-deps})
set_par_sizes(test-mpi-cgsolve-2d 4)

configure_file(test-cgsolve-0.json ${CMAKE_BINARY_DIR}/test-cgsolve-0.json COPYONLY)
configure_file(test-cgsolve-1.json ${CMAKE_BINARY_DIR}/test-cgsolve-1.json COPYONLY)

add_par_unit(test-mpi-cgredist-2d ${CMAKE_CURRENT_SOURCE_DIR}/test_cgredist.cc ${CMAKE_CURRENT_SOURCE_DIR}/../../mpi_main.cc)
target_link_libraries(test-mpi-cgredist-2d cedar ${cedar-deps})
set_par_sizes(test-mpi-cgredist-2d 16)

configure_file(test-cgredist-0.json ${CMAKE_BINARY_DIR}/test-cgredist-0.json COPYONLY)
configure_file(test-cgredist-1.json ${CMAKE_BINARY_DIR}/test-cgredist-1.json COPYONLY)
