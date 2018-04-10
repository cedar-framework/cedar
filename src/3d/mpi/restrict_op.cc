#include <cedar/3d/mpi/prolong_op.h>
#include <cedar/3d/mpi/restrict_op.h>
#include <cedar/3d/mpi/grid_func.h>


namespace cedar { namespace cdr3 { namespace mpi {

std::ostream & operator<< (std::ostream &os, const restrict_op &R)
{
	os << *R.P;
	return os;
}

}}}
