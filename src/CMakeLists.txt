set(dirs 2d 3d config mpi perf util)

foreach(dir ${dirs})
  add_subdirectory(${dir})
  foreach(dsrc ${${dir}_SOURCES})
	list(APPEND sources "src/${dir}/${dsrc}")
  endforeach(dsrc)
endforeach(dir)

set(boxmg_src ${sources} PARENT_SCOPE)
