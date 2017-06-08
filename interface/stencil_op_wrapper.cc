#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <cedar/2d/stencil_op.h>



namespace py = pybind11;
using namespace cedar;
using namespace cedar::cdr2;

template<class sten>
void export_stencil_op(py::module & m, std::string sten_name) {
	std::string name("stencil_op_" + sten_name);
	py::class_<stencil_op<sten>>(m, name.c_str(), py::buffer_protocol())
		.def(py::init<len_t,len_t>())
		.def("__repr__", [](stencil_op<sten> & so) -> std::string {
				std::stringstream ss;
				ss << so;
				return ss.str();
			})
		.def_buffer([](stencil_op<sten> &so) -> py::buffer_info {
				return py::buffer_info(
					so.data(),
					sizeof(real_t),
					py::format_descriptor<real_t>::format(),
					3,
					{so.len(0), so.len(1), stencil_ndirs<sten>::value},
					{sizeof(real_t),
							sizeof(real_t)*so.len(0),
							sizeof(real_t)*so.len(0)*so.len(1)}
					);
			});
}

PYBIND11_PLUGIN(pycedar) {
	py::module m("pycedar", "");
	export_stencil_op<five_pt>(m,"five");
	export_stencil_op<nine_pt>(m,"nine");


	return m.ptr();
}

