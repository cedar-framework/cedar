pybind11_add_module(_cdr2 wrap2.cc
  wrap2_gridfunc.cc)

target_link_libraries(_cdr2 cedar ${cedar-deps})


pybind11_add_module(gallery wrap2_gallery.cc)

target_link_libraries(gallery cedar ${cedar-deps})
