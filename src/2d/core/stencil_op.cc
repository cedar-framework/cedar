#include <iomanip>

#include "kernel/name.h"
#include "stencil_op.h"

using namespace boxmg::bmg2d;


void StencilOp::set_registry(std::shared_ptr<KernelRegistry> kreg)
{
	kernels = kreg;
}


std::shared_ptr<boxmg::KernelRegistry> StencilOp::get_registry() const
{
	return kernels;
}


void StencilOp::apply(const GridFunc &x, GridFunc &b) const
{
	if (kernels != nullptr) {
		kernels->run(kernel::name::matvec,
		             static_cast<const StencilOp&>(*this),
		             x,b);
	} else {
		log::error << "Kernel registry for stencil op not set!" << std::endl;
	}

}


void StencilOp::residual(const GridFunc &x, const GridFunc &b, GridFunc & r) const
{
	if (kernels != nullptr) {
		kernels->run(kernel::name::residual,
		             static_cast<const StencilOp&>(*this),
		             x,b,r);
	} else {
		log::error << "Kernel registry for stencil op not set!" << std::endl;
	}
}


GridFunc StencilOp::residual(const GridFunc &x, const GridFunc &b) const
{
	auto r = GridFunc(x.shape(0),x.shape(1));
	this->residual(x,b,r);

	return r;
}


namespace boxmg { namespace bmg2d {

std::ostream & operator<< (std::ostream & os, const StencilOp &op)
{
	auto sten = op.stencil();
	unsigned int width = 4;

	os << std::setprecision(7);

	if (sten.five_pt()) {
		for (auto j: sten.range(1)) {
			for (auto i: sten.range(0)) {
				os << std::setw(width) << (j+1-2)*sten.len(0) + i+1-2 << " "
				   << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
				   << std::scientific << -sten(i,j,Dir::N)
				   << std::scientific << -sten(i,j,Dir::W) << " "
				   << std::scientific <<  sten(i,j,Dir::C)
				   << std::scientific << -sten(i,j,Dir::E)
				   << std::scientific << -sten(i,j,Dir::S) << '\n';
			}
		}
	} else {
		for (auto j: sten.range(1)) {
			for (auto i: sten.range(0)) {
				os << std::setw(width) << (j+1-2)*sten.len(0) + i+1-2 << " "
				   << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
				   << std::scientific << -sten(i,j,Dir::NW)
				   << std::scientific << -sten(i,j,Dir::N)
				   << std::scientific << -sten(i,j,Dir::NE)
				   << std::scientific << -sten(i,j,Dir::W) << " "
				   << std::scientific <<  sten(i,j,Dir::C)
				   << std::scientific << -sten(i,j,Dir::E)
				   << std::scientific << -sten(i,j,Dir::SW)
				   << std::scientific << -sten(i,j,Dir::S)
				   << std::scientific << -sten(i,j,Dir::SE) << '\n';
			}
		}
	}

	return os;
}

}}
