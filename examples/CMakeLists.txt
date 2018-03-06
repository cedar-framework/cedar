add_executable(ser-poisson-2d basic-2d-ser/poisson.cc)
target_link_libraries(ser-poisson-2d cedar ${cedar-deps})

add_executable(ser-periodic-2d basic-2d-ser/periodic.cc)
target_link_libraries(ser-periodic-2d cedar ${cedar-deps})

if (ENABLE_MPI)
  add_executable(mpi-poisson-2d basic-2d-mpi/poisson.cc)
  target_link_libraries(mpi-poisson-2d cedar ${cedar-deps})

  add_executable(mpi-periodic-2d basic-2d-mpi/periodic.cc)
  target_link_libraries(mpi-periodic-2d cedar ${cedar-deps})

  add_executable(periodic-halo-2d halo-2d/per_halo.cc)
  target_link_libraries(periodic-halo-2d cedar ${cedar-deps})

  if (ENABLE_3D)
	add_executable(mpi-poisson-3d basic-3d-mpi/poisson.cc)
	target_link_libraries(mpi-poisson-3d cedar ${cedar-deps})

	add_executable(mpi-periodic-3d basic-3d-mpi/periodic.cc)
	target_link_libraries(mpi-periodic-3d cedar ${cedar-deps})

	add_executable(periodic-halo-3d halo-3d/per_halo.cc)
	target_link_libraries(periodic-halo-3d cedar ${cedar-deps})
  endif()
endif()

if (ENABLE_3D)
  add_executable(ser-poisson-3d basic-3d-ser/poisson.cc)
  target_link_libraries(ser-poisson-3d cedar ${cedar-deps})

  add_executable(ser-periodic-3d basic-3d-ser/periodic.cc)
  target_link_libraries(ser-periodic-3d cedar ${cedar-deps})
endif()

configure_file(basic-2d-mpi/poisson-config.json ${CMAKE_BINARY_DIR}/config.json COPYONLY)
