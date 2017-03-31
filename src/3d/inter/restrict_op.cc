#include <cedar/kernel_name.h>
#include <cedar/3d/inter/restrict_op.h>

namespace cedar { namespace cdr3 { namespace inter {

std::ostream & operator<< (std::ostream &os, const restrict_op &R)
{
	os << *R.P;
	return os;
}


void restrict_op::apply(const grid_func &x, grid_func &y) const
{
	if (kernels != nullptr) {
		kernels->run(kernel_name::restriction,
		             static_cast<const restrict_op&>(*this),
		             x,y);
	} else {
		log::error << "You forgot to give the Restrict Operator a kernel registry!" << std::endl;
	}
}


grid_func operator*(const restrict_op & R, const grid_func &x)
{
	auto & P = R.getP();
	auto & Psten = P.stencil();

	grid_func y(Psten.shape(0), Psten.shape(1), Psten.shape(2));
	R.apply(x,y);

	return y;
}

}}}
