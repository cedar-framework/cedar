if (ENABLE_UNIT_TESTS)
  enable_testing()
  find_package(GTest REQUIRED)
  if (GTEST_FOUND)
	include_directories(${GTEST_INCLUDE_DIRS})
  else()
	add_subdirectory(${CMAKE_SOURCE_DIR}/external/gtest-1.7.0)
	include_directories(${CMAKE_SOURCE_DIR}/external/gtest-1.7.0/include)
	set(GTEST_BOTH_LIBRARIES gtest gtest_main)
  endif()
endif()


function(add_par_unit target sources)
  set(src ${ARGV})
  list(REMOVE_AT src 0)
  add_executable(${target} ${src})
  target_link_libraries(${target} ${GTEST_BOTH_LIBRARIES} ${boxmg-deps} boxmg)
  set_target_properties(${target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test)
endfunction(add_par_unit)


function(set_par_sizes target sizes)
  set(comm_sizes ${ARGV})
  list(REMOVE_AT comm_sizes 0)
  foreach(commsize ${comm_sizes})
	set(test_parameters -np ${commsize} "${CMAKE_BINARY_DIR}/test/${target}")
	add_test(NAME ${target}-${commsize} COMMAND "mpiexec" ${test_parameters})
  endforeach()
endfunction()


function(add_unit target sources)
  set(src ${ARGV})
  list(REMOVE_AT src 0)
  add_executable(${target} ${src})
  target_link_libraries(${target} ${GTEST_BOTH_LIBRARIES} ${boxmg-deps} boxmg)
  set_target_properties(${target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test)
  if(ENABLE_JENKINS_OUTPUT)
    add_test(${target} ${CMAKE_BINARY_DIR}/test/${target}
      --gtest_output=xml:${CMAKE_BINARY_DIR}/test/${target}.xml
      --gtest_color=yes)
  else()
	add_test(NAME ${target} COMMAND ${CMAKE_BINARY_DIR}/test/${target})
  endif()
endfunction(add_unit)
