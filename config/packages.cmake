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

if(ENABLE_MPI)
	set(boxmg-examples_LINKER_FLAGS ${Boost_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${MPI_Fortran_LIBRARIES})
endif(ENABLE_MPI)


include_directories(${CMAKE_BINARY_DIR}/include/boxmg/2d)
