add_executable(solve-test solve/solve_test.cc)
target_link_libraries(solve-test boxmg ${boxmg-deps})

add_executable(mpi-test mpi/solve_test.cc mpi/diffusion.f90 mpi/putf.f90)
target_link_libraries(mpi-test boxmg ${boxmg-deps})

add_executable(solve-test-3 solve-3d/test.cc solve-3d/putf.f90 solve-3d/put_sol.f90)
target_link_libraries(solve-test-3 boxmg ${boxmg-deps})

add_executable(mpi-test-3 mpi-3d/solve_test.cc mpi-3d/putf.f90 mpi-3d/put_sol.f90)
target_link_libraries(mpi-test-3 boxmg ${boxmg-deps})

# add_executable(timer-test timer-test/timer_test.cc)
# target_link_libraries(timer-test boxmg ${boxmg-deps})
