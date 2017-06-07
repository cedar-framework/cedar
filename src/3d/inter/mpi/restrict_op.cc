#include <cedar/kernel_name.h>
#include <cedar/3d/inter/mpi/prolong_op.h>
#include <cedar/3d/inter/mpi/restrict_op.h>
#include <cedar/3d/mpi/grid_func.h>


namespace cedar { namespace cdr3 { namespace inter { namespace mpi {

std::ostream & operator<< (std::ostream &os, const restrict_op &R)
{
	os << *R.P;
	return os;
}

}}}}
