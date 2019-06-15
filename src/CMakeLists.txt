set(dirs 2d perf util interface services)

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
  src/global_params.cc
  src/global_manager.cc
  src/config.cc
  src/multilevel_settings.cc
  src/malloc_pool.cc
  src/memory.cc
  src/log.cc
  src/init.cc
  PARENT_SCOPE)
