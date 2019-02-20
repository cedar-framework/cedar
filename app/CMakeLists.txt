add_executable(openmp-relax relax.cc)
target_link_libraries(openmp-relax cedar)
configure_file(openmp-relax.json ${CMAKE_BINARY_DIR}/openmp-relax.json)
