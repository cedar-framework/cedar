#include <iomanip>

#include <cedar/2d/stencil_op.h>

using namespace cedar::cdr2;

namespace cedar { namespace cdr2 {

std::ostream & operator<< (std::ostream & os, const stencil_op &op)
{
	auto & sten = op.stencil();
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
