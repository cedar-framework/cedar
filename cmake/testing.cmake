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
  set(mpimain_SRC ${CMAKE_SOURCE_DIR}/test/mpi_main.cc)
  add_executable(${target} ${sources} ${mpimain_SRC})
  target_link_libraries(${target} ${GTEST_BOTH_LIBRARIES} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} efield-test)
  set_target_properties(${target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test)
  if(ENABLE_JENKINS_OUTPUT)
    add_test(${target} ${CMAKE_BINARY_DIR}/test/${target}
      --gtest_output=xml:${CMAKE_BINARY_DIR}/test/${target}.xml
      --gtest_color=yes)
  else()
	if (${ARGC} GREATER 2)
	  set(commsize ${ARGV2})
	else()
	  set(commsize ${TEST_COMM_SIZE})
	endif()
	set(test_parameters -np ${commsize} "${CMAKE_BINARY_DIR}/test/${target}")
	add_test(NAME ${target} COMMAND "mpiexec" ${test_parameters})
  endif()
endfunction(add_par_unit)


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
