#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <cedar/2d/gallery.h>

namespace py = pybind11;
using namespace cedar;
using namespace cedar::cdr2;


PYBIND11_PLUGIN(gallery) {
	py::module m("gallery", "");

	m.def("poisson", &gallery::poisson);

	return m.ptr();
}
