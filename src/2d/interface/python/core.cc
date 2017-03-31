#include <boost/python.hpp>
#include "core/grid_stencil.h"
#include "core/grid_func.h"
#include "core/types.h"


BOOST_PYTHON_MODULE(core)
{
	using namespace boost::python;
	using namespace cedar;
	using namespace cedar::cdr2::core;
	class_<GridStencil>("GridStencil", init<len_t,len_t,unsigned int, bool, bool, bool>());
	enum_<Dir>("Dir")
		.value("C", Dir::C)
		.value("N", Dir::N)
		.value("S", Dir::S)
		.value("E", Dir::E)
		.value("W", Dir::W);
	//class_<GridFunc>("GridFunc", init<len_t, len_t>());
}
