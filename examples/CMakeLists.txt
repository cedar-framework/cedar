add_executable(solve-test solve/solve_test.cc)
target_link_libraries(solve-test cedar ${cedar-deps})

if (ENABLE_MPI)
  add_executable(mpi-test mpi/solve_test.cc)
  target_link_libraries(mpi-test cedar ${cedar-deps})

  add_executable(halo-2d halo/halo_test.cc)
  target_link_libraries(halo-2d cedar ${cedar-deps})

  if (ENABLE_3D)
	add_executable(mpi-test-3 mpi-3d/solve_test.cc)
	target_link_libraries(mpi-test-3 cedar ${cedar-deps})
  endif()
endif()

if (ENABLE_3D)
  add_executable(solve-test-3 solve-3d/test.cc)
  target_link_libraries(solve-test-3 cedar ${cedar-deps})
endif()

# add_executable(timer-test timer-test/timer_test.cc)
# target_link_libraries(timer-test cedar ${cedar-deps})

configure_file(mpi/config.json ${CMAKE_BINARY_DIR}/config.json COPYONLY)
configure_file(halo/config.json ${CMAKE_BINARY_DIR}/halo-config.json COPYONLY)
