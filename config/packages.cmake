include(cxx11)

check_for_cxx11_compiler(CXX11_COMPILER)

if(CXX11_COMPILER)
  enable_cxx11()
else()
  message(FATAL_ERROR "C++11 compatible compiler not found")
endif()

FIND_PACKAGE(Boost 1.44 REQUIRED COMPONENTS system filesystem)
include_directories(${Boost_INCLUDE_DIR})

FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

find_package(MPI REQUIRED)
include_directories(${MPI_C_INCLUDE_PATH})

set(boxmg-examples_LINKER_FLAGS ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES})

include_directories(${CMAKE_BINARY_DIR}/include/boxmg/2d)
