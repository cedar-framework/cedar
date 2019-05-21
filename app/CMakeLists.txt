add_executable(restest restest.cc)
target_link_libraries(restest cedar)
configure_file(restest.json ${CMAKE_BINARY_DIR}/restest.json COPYONLY)
