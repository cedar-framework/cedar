# add_executable(test-thread test.cc)
# target_compile_options(test-thread PUBLIC -fopenmp -I/home/andrew/local/include -L/home/andrew/local/lib -lomp -Wl,-rpath=/home/andrew/local/lib)
# target_link_libraries(test-thread cedar -fopenmp -I/home/andrew/local/include -L/home/andrew/local/lib -lomp -Wl,-rpath=/home/andrew/local/lib)
# target_compile_options(test-thread PUBLIC -I/home/andrew/local/include)
# target_link_libraries(test-thread cedar -L/home/andrew/local/lib -labt)


# add_executable(test-coll coll.cc)
# target_link_libraries(test-coll cedar)


# add_executable(plane-run plane.cc)
# target_link_libraries(plane-run cedar)

add_executable(line-run line.cc line.f90)
target_link_libraries(line-run cedar)
