add_par_unit(test-mpi-poisson-3d ${CMAKE_CURRENT_SOURCE_DIR}/test_poisson.cc ${CMAKE_CURRENT_SOURCE_DIR}/../../mpi_main.cc)
target_link_libraries(test-mpi-poisson-3d cedar)
set_par_sizes(test-mpi-poisson-3d 1 3 4 8 9)


add_par_unit(test-mpi-relax-3d ${CMAKE_CURRENT_SOURCE_DIR}/test_relax.cc ${CMAKE_CURRENT_SOURCE_DIR}/../../mpi_main.cc)
target_link_libraries(test-mpi-relax-3d cedar)
set_par_sizes(test-mpi-relax-3d 1 3 8 9)


add_par_unit(test-mpi-planes-3d ${CMAKE_CURRENT_SOURCE_DIR}/test_planes.cc ${CMAKE_CURRENT_SOURCE_DIR}/../../mpi_main.cc)
target_link_libraries(test-mpi-planes-3d cedar)
set_par_sizes(test-mpi-planes-3d 1 3)# 8 9)
configure_file(test-planes-mpi.json ${CMAKE_BINARY_DIR}/test-planes-mpi.json COPYONLY)


add_par_unit(test-mpi-cgredist-3d ${CMAKE_CURRENT_SOURCE_DIR}/test_cgredist.cc ${CMAKE_CURRENT_SOURCE_DIR}/../../mpi_main.cc)
target_link_libraries(test-mpi-cgredist-3d cedar)
set_par_sizes(test-mpi-cgredist-3d 8)

if (ENABLE_PLANE_AGG)
  add_par_unit(test-mpi-planeagg ${CMAKE_CURRENT_SOURCE_DIR}/test_planeagg.cc ${CMAKE_CURRENT_SOURCE_DIR}/../../abt_main.cc)
  target_link_libraries(test-mpi-planeagg cedar)
  set_par_sizes(test-mpi-planeagg 3 4 8 16)
endif()

configure_file(test-cgredist-0.json ${CMAKE_BINARY_DIR}/test-cgredist3-0.json COPYONLY)
configure_file(test-cgredist-1.json ${CMAKE_BINARY_DIR}/test-cgredist3-1.json COPYONLY)
configure_file(test-planeagg-lines.json ${CMAKE_BINARY_DIR}/test-planeagg-lines.json COPYONLY)
configure_file(test-planeagg.json ${CMAKE_BINARY_DIR}/test-planeagg.json COPYONLY)
