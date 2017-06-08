pybind11_add_module(testpy test.cc)
pybind11_add_module(pycedar stencil_op_wrapper.cc)
target_link_libraries(pycedar cedar ${cedar-deps})
