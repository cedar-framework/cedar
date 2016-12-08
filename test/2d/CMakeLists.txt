add_unit(residual-2d ${CMAKE_CURRENT_SOURCE_DIR}/test_residual.cc)
target_link_libraries(residual-2d boxmg ${boxmg-deps})
add_unit(relax-2d ${CMAKE_CURRENT_SOURCE_DIR}/pyrelax.c ${CMAKE_CURRENT_SOURCE_DIR}/test_relax.cc)
target_link_libraries(relax-2d boxmg ${boxmg-deps})
