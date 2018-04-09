#include <cedar/kernel_name.h>
#include <cedar/2d/mpi/prolong_op.h>
#include <cedar/2d/mpi/restrict_op.h>
#include <cedar/2d/mpi/grid_func.h>


namespace cedar { namespace cdr2 { namespace mpi {

std::ostream & operator<<(std::ostream & os, const restrict_op & R)
{
	os << *R.P;
	return os;
}

}}}
