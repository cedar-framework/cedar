#include <cedar/kernel_name.h>

#include <cedar/2d/inter/restrict_op.h>


namespace cedar { namespace cdr2 { namespace inter {

std::ostream & operator<< (std::ostream &os, const restrict_op &R)
{
	os << *R.P;
	return os;
}

}}}
