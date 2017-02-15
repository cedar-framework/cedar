cmake_minimum_required(VERSION 2.8)
project(boxmg)

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

set(boxmg-deps ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES})

include_directories(${CMAKE_BINARY_DIR}/include/boxmg/2d)

add_subdirectory(src)
add_library(boxmg ${boxmg_src})
install(TARGETS boxmg LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

if(ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

if(ENABLE_UNIT_TESTS)
  find_package(PythonLibs REQUIRED)
  find_package(NumPy REQUIRED)
  message(STATUS ${PYTHON_INCLUDE_DIR})
  list(APPEND boxmg-deps ${PYTHON_LIBRARY})
  include_directories(${PYTHON_INCLUDE_DIR} ${PYTHON_NUMPY_INCLUDE_DIR})

  add_subdirectory(test)
endif()
