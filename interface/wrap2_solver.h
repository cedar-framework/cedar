#ifndef WRAP2_SOLVER_H
#define WRAP2_SOLVER_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <cedar/2d/solver.h>

namespace py = pybind11;
using namespace cedar;
using namespace cedar::cdr2;

template<class sten>
void export_solver2(py::module & m, std::string sten_name)
{
	m.def("solve", &solve<five_pt>);
	m.def("solve", &solve<nine_pt>);
}


#endif
