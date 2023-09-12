#include <cedar/2d/gpu/prolong_op.h>
#include <cedar/2d/gpu/restrict_op.h>
#include <cedar/2d/gpu/grid_func.h>


namespace cedar { namespace cdr2 { namespace gpu {namespace mpi {

std::ostream & operator<<(std::ostream & os, const restrict_op & R)
{
	os << *R.P;
	return os;
}

}}}}
