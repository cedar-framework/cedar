#include <cedar/kernel_name.h>
#include <cedar/2d/inter/mpi/prolong_op.h>
#include <cedar/2d/inter/mpi/restrict_op.h>
#include <cedar/2d/mpi/grid_func.h>


namespace cedar { namespace cdr2 { namespace inter { namespace mpi {
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
