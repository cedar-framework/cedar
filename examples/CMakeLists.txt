add_executable(solve-test solve/solve_test.cc)
target_link_libraries(solve-test boxmg ${boxmg-deps})

if (ENABLE_MPI)
  add_executable(mpi-test mpi/solve_test.cc)
  target_link_libraries(mpi-test boxmg ${boxmg-deps})

  add_executable(halo-2d halo/halo_test.cc)
  target_link_libraries(halo-2d boxmg ${boxmg-deps})

  if (ENABLE_3D)
	add_executable(mpi-test-3 mpi-3d/solve_test.cc)
	target_link_libraries(mpi-test-3 boxmg ${boxmg-deps})
  endif()
endif()

if (ENABLE_3D)
  add_executable(solve-test-3 solve-3d/test.cc)
  target_link_libraries(solve-test-3 boxmg ${boxmg-deps})
endif()

# add_executable(timer-test timer-test/timer_test.cc)
# target_link_libraries(timer-test boxmg ${boxmg-deps})

configure_file(mpi/config.json ${CMAKE_BINARY_DIR}/config.json COPYONLY)
configure_file(halo/config.json ${CMAKE_BINARY_DIR}/halo-config.json COPYONLY)
