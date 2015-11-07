#include "kernel/name.h"

#include "restrict_op.h"


namespace boxmg { namespace bmg2d { namespace inter {

std::ostream & operator<< (std::ostream &os, const RestrictOp &R)
{
	os << *R.P;
	return os;
}



void RestrictOp::set_registry(std::shared_ptr<KernelRegistry> kreg)
{
	kernels = kreg;
}


void RestrictOp::apply(const GridFunc &x, GridFunc &y) const
{
	if (kernels != nullptr) {
		kernels->run(kernel::name::restriction,
		             static_cast<const RestrictOp&>(*this),
		             x,y);
	} else {
		log::error << "You forgot to give the Restrict Operator a kernel registry!" << std::endl;
	}
}


GridFunc operator*(const RestrictOp & R, const GridFunc &x)
{
	auto & P = R.getP();
	auto & Psten = P.stencil();

	GridFunc y(Psten.shape(0), Psten.shape(1));
	R.apply(x,y);

	return y;
}

}}}
