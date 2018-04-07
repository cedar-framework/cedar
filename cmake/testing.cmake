if (ENABLE_UNIT_TESTS)
  enable_testing()
  find_package(GTest QUIET)
  if (GTEST_FOUND)
	include_directories(${GTEST_INCLUDE_DIRS})
  else()
	# Download and unpack googletest at configure time
	configure_file(cmake/testing.in googletest-download/CMakeLists.txt)
	execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
	  RESULT_VARIABLE result
	  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
	if(result)
	  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
	endif()
	execute_process(COMMAND ${CMAKE_COMMAND} --build .
	  RESULT_VARIABLE result
	  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
	if(result)
	  message(FATAL_ERROR "Build step for googletest failed: ${result}")
	endif()

	# Prevent overriding the parent project's compiler/linker
	# settings on Windows
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

	# Add googletest directly to our build. This defines
	# the gtest and gtest_main targets.
	add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
      ${CMAKE_BINARY_DIR}/googletest-build
      EXCLUDE_FROM_ALL)

	# The gtest/gtest_main targets carry header search path
	# dependencies automatically when using CMake 2.8.11 or
	# later. Otherwise we have to add them here ourselves.
	if (CMAKE_VERSION VERSION_LESS 2.8.11)
	  include_directories("${gtest_SOURCE_DIR}/include")
	endif()

	set(GTEST_BOTH_LIBRARIES gtest gtest_main)
  endif()
endif()


function(add_par_unit target sources)
  set(src ${ARGV})
  list(REMOVE_AT src 0)
  add_executable(${target} ${src})
  target_link_libraries(${target} ${GTEST_BOTH_LIBRARIES} ${cedar-deps} cedar)
  set_target_properties(${target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test)
endfunction(add_par_unit)


function(set_par_sizes target sizes)
  set(comm_sizes ${ARGV})
  list(REMOVE_AT comm_sizes 0)
  foreach(commsize ${comm_sizes})
	set(test_parameters -np ${commsize} "${CMAKE_BINARY_DIR}/test/${target}")
	add_test(NAME ${target}-${commsize} COMMAND "mpiexec" ${test_parameters} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  endforeach()
endfunction()


function(add_unit target sources)
  set(src ${ARGV})
  list(REMOVE_AT src 0)
  add_executable(${target} ${src})
  target_link_libraries(${target} ${GTEST_BOTH_LIBRARIES} ${cedar-deps} cedar)
  set_target_properties(${target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test)
  if(ENABLE_JENKINS_OUTPUT)
    add_test(${target} ${CMAKE_BINARY_DIR}/test/${target}
      --gtest_output=xml:${CMAKE_BINARY_DIR}/test/${target}.xml
      --gtest_color=yes)
  else()
	add_test(NAME ${target} COMMAND ${CMAKE_BINARY_DIR}/test/${target} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  endif()
endfunction(add_unit)
