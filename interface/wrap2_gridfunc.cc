#include "wrap2_gridfunc.h"


void export_gridfunc2(py::module & m) {
	py::class_<grid_func>(m, "grid_func", py::buffer_protocol())
		.def(py::init<len_t,len_t>())
		.def("__repr__", [](grid_func & f) -> std::string {
				std::stringstream ss;
				ss << f;
				return ss.str();
			})
		.def_buffer([](grid_func & f) -> py::buffer_info {
				return py::buffer_info(
					f.data(),
					sizeof(real_t),
					py::format_descriptor<real_t>::format(),
					2,
					{f.len(0), f.len(1)},
					{sizeof(real_t), sizeof(real_t)*f.len(0)}
					);
			});
}

