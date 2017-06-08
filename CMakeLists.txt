cmake_minimum_required(VERSION 2.8)
project(cedar)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(cxx11)
include(cxx_config)
include(options)
include(testing)

enable_language(C CXX Fortran)

add_definitions(-DMSG_DOUBLE -DNO_SECOND_UNDERSCORE -fPIC -DBOUNDS_CHECK)
include_directories(include)
install(DIRECTORY include/ DESTINATION include)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)


FIND_PACKAGE(Boost 1.44 REQUIRED COMPONENTS system filesystem)
include_directories(${Boost_INCLUDE_DIR})

FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

find_package(MPI REQUIRED)
include_directories(${MPI_C_INCLUDE_PATH})
include_directories(${MPI_Fortran_INCLUDE_PATH})

set(cedar-deps ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES})

add_subdirectory(src)
add_library(cedar ${cedar_src})
install(TARGETS cedar LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

if(ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

if(ENABLE_UNIT_TESTS)
  find_package(PythonLibs 2.7 REQUIRED)
  find_package(NumPy REQUIRED)
  message(STATUS ${PYTHON_INCLUDE_DIR})
  list(APPEND cedar-deps ${PYTHON_LIBRARY})
  include_directories(${PYTHON_INCLUDE_DIR} ${PYTHON_NUMPY_INCLUDE_DIR})

  add_subdirectory(test)
endif()

add_subdirectory(pybind11)
add_subdirectory(interface)
