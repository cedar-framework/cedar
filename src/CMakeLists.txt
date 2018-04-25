set(dirs 2d perf util interface)

if (ENABLE_3D)
  list(APPEND dirs 3d)
endif()

if(ENABLE_MPI)
  list(APPEND dirs mpi)
endif()

foreach(dir ${dirs})
  add_subdirectory(${dir})
  foreach(dsrc ${${dir}_SOURCES})
	list(APPEND sources "src/${dir}/${dsrc}")
  endforeach(dsrc)
endforeach(dir)

set(cedar_src ${sources}
  src/kernel_params.cc
  src/config.cc
  src/multilevel_settings.cc
  PARENT_SCOPE)
