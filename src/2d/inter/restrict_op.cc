#include <boxmg/2d/kernel/name.h>

#include <boxmg/2d/inter/restrict_op.h>


namespace boxmg { namespace bmg2d { namespace inter {

std::ostream & operator<< (std::ostream &os, const restrict_op &R)
{
	os << *R.P;
	return os;
}



void restrict_op::set_registry(std::shared_ptr<KernelRegistry> kreg)
{
	kernels = kreg;
}


void restrict_op::apply(const grid_func &x, grid_func &y) const
{
	if (kernels != nullptr) {
		kernels->run(kernel::name::restriction,
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

	grid_func y(Psten.shape(0), Psten.shape(1));
	R.apply(x,y);

	return y;
}

}}}
