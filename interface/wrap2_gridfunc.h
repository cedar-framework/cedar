#ifndef WRAP2_GRIDFUNC_H
#define WRAP2_GRIDFUNC_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <cedar/2d/grid_func.h>

namespace py = pybind11;
using namespace cedar;
using namespace cedar::cdr2;

void export_gridfunc2(py::module & m);

#endif
