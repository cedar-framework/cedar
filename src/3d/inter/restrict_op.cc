#include <cedar/kernel_name.h>
#include <cedar/3d/inter/restrict_op.h>

namespace cedar { namespace cdr3 { namespace inter {

std::ostream & operator<< (std::ostream &os, const restrict_op &R)
{
	os << *R.P;
	return os;
}

}}}
