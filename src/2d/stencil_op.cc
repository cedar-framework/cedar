#include <iomanip>

#include <boxmg/2d/kernel/name.h>
#include <boxmg/2d/stencil_op.h>

using namespace boxmg::bmg2d;


void stencil_op::set_registry(std::shared_ptr<KernelRegistry> kreg)
{
	kernels = kreg;
}


std::shared_ptr<boxmg::KernelRegistry> stencil_op::get_registry() const
{
	return kernels;
}


void stencil_op::apply(const grid_func &x, grid_func &b) const
{
	if (kernels != nullptr) {
		kernels->run(kernel::name::matvec,
		             static_cast<const stencil_op&>(*this),
		             x,b);
	} else {
		log::error << "Kernel registry for stencil op not set!" << std::endl;
	}

}


void stencil_op::residual(const grid_func &x, const grid_func &b, grid_func & r) const
{
	if (kernels != nullptr) {
		kernels->run(kernel::name::residual,
		             static_cast<const stencil_op&>(*this),
		             x,b,r);
	} else {
		log::error << "Kernel registry for stencil op not set!" << std::endl;
	}
}


grid_func stencil_op::residual(const grid_func &x, const grid_func &b) const
{
	auto r = grid_func(x.shape(0),x.shape(1));
	this->residual(x,b,r);

	return r;
}


namespace boxmg { namespace bmg2d {

std::ostream & operator<< (std::ostream & os, const stencil_op &op)
{
	auto sten = op.stencil();
	unsigned int width = 4;

	os << std::setprecision(7);

	if (sten.five_pt()) {
		for (auto j: sten.range(1)) {
			for (auto i: sten.range(0)) {
				os << std::setw(width) << (j+1-2)*sten.len(0) + i+1-2 << " "
				   << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
				   << std::scientific << -sten(i,j,dir::N)
				   << std::scientific << -sten(i,j,dir::W) << " "
				   << std::scientific <<  sten(i,j,dir::C)
				   << std::scientific << -sten(i,j,dir::E)
				   << std::scientific << -sten(i,j,dir::S) << '\n';
			}
		}
	} else {
		for (auto j: sten.range(1)) {
			for (auto i: sten.range(0)) {
				os << std::setw(width) << (j+1-2)*sten.len(0) + i+1-2 << " "
				   << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
				   << std::scientific << -sten(i,j,dir::NW)
				   << std::scientific << -sten(i,j,dir::N)
				   << std::scientific << -sten(i,j,dir::NE)
				   << std::scientific << -sten(i,j,dir::W) << " "
				   << std::scientific <<  sten(i,j,dir::C)
				   << std::scientific << -sten(i,j,dir::E)
				   << std::scientific << -sten(i,j,dir::SW)
				   << std::scientific << -sten(i,j,dir::S)
				   << std::scientific << -sten(i,j,dir::SE) << '\n';
			}
		}
	}

	return os;
}

}}
