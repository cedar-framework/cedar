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

mpi::grid_func operator*(const restrict_op & R, const mpi::grid_func &x)
{
	auto & P = R.getP();
	mpi::grid_func y(P.grid_ptr());

	R.apply(x,y);

	return y;
}


void restrict_op::apply(const mpi::grid_func &x, mpi::grid_func &y) const
{
	if (kernels != nullptr) {
		kernels->run(kernel_name::restriction,
		             static_cast<const restrict_op&>(*this),
		             x,y);
	} else {
		log::error << "You forgot to give the Restrict Operator a kernel registry!" << std::endl;
	}
}

}}}}
