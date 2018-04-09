#include <cedar/2d/restrict_op.h>


namespace cedar { namespace cdr2 {

std::ostream & operator<< (std::ostream &os, const restrict_op &R)
{
	os << *R.P;
	return os;
}

}}
