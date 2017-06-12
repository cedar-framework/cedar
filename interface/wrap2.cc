#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "wrap2_stencil.h"
#include "wrap2_gridfunc.h"
#include "wrap2_solver.h"

PYBIND11_PLUGIN(_cdr2) {
	py::module m("_cdr2", "");
	export_stencil_op2<five_pt>(m, "five");
	export_stencil_op2<nine_pt>(m, "nine");
	export_solver2<five_pt>(m, "five");
	export_solver2<nine_pt>(m, "nine");
	export_gridfunc2(m);

	return m.ptr();
}
