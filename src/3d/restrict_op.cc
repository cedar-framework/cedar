#include <cedar/3d/restrict_op.h>

namespace cedar { namespace cdr3 {

std::ostream & operator<< (std::ostream &os, const restrict_op &R)
{
	os << *R.P;
	return os;
}

}}
