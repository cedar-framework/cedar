set(dirs 2d config perf util)

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

set(boxmg_src ${sources} PARENT_SCOPE)
