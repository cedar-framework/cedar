#include <iomanip>
#include <cedar/3d/stencil_op.h>


namespace cedar { namespace cdr3 {

std::ostream & operator<< (std::ostream & os, const stencil_op & op)
{
	auto & sten = op.stencil();
	auto nx = sten.len(0);
	auto ny = sten.len(1);

	unsigned int width = 4;

	os << std::setprecision(7);
	if (sten.five_pt()) {
		for (auto k : sten.range(2)) {
			for (auto j : sten.range(1)) {
				for (auto i : sten.range(0)) {
					os << std::setw(width) << (k+1-2)*(nx*ny) + (j+1-2)*nx + i+1-2 << " "
					   << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
					   << std::setw(width) << k + 1 << ", "
					   << std::scientific << -sten(i,j+1,k,dir::PS) << ", "
					   << std::scientific << -sten(i,j,k+1,dir::B) << ", "
					   << std::scientific << -sten(i,j,k,dir::PW) << ",  "
					   << std::scientific << sten(i,j,k,dir::P) << ", "
					   << std::scientific << -sten(i+1,j,k,dir::PW) << ", "
					   << std::scientific << -sten(i,j,k,dir::B) << ", "
					   << std::scientific << -sten(i,j,k,dir::PS) << '\n';
				}
			}
		}
	} else {
		for (auto k : sten.range(2)) {
			for (auto j : sten.range(1)) {
				for (auto i : sten.range(0)) {
					os << std::setw(width) << (k+1-2)*(nx*ny) + (j+1-2)*nx + i+1-2 << " "
					   << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
					   << std::setw(width) << k + 1 << ", " << " N, "
					   << std::scientific << -sten(i,j+1,k+1,dir::BSE) << ", "
					   << std::scientific << -sten(i,j+1,k+1,dir::BS) << ", "
					   << std::scientific << -sten(i+1,j+1,k+1,dir::BSW) << ", "
					   << std::scientific << -sten(i,j+1,k,dir::PNW) << ", "
					   << std::scientific << -sten(i,j+1,k,dir::PS) << ", "
					   << std::scientific << -sten(i+1,j+1,k,dir::PSW) << ", "
					   << std::scientific << -sten(i,j+1,k,dir::BNW) << ", "
					   << std::scientific << -sten(i,j+1,k,dir::BN) << ", "
					   << std::scientific << -sten(i+1,j+1,k,dir::BNE) << '\n';

					os << std::setw(width) << (k+1-2)*(nx*ny) + (j+1-2)*nx + i+1-2 << " "
					   << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
					   << std::setw(width) << k + 1 << ", " << " O, "
					   << std::scientific << -sten(i,j,k+1,dir::BE) << ", "
					   << std::scientific << -sten(i,j,k+1,dir::B) << ", "
					   << std::scientific << -sten(i+1,j,k+1,dir::BW) << ", "
					   << std::scientific << -sten(i,j,k,dir::PW) << ", "
					   << std::scientific << sten(i,j,k,dir::P) << ", "
					   << std::scientific << -sten(i+1,j,k,dir::PW) << ", "
					   << std::scientific << -sten(i,j,k,dir::BW) << ", "
					   << std::scientific << -sten(i,j,k,dir::B) << ", "
					   << std::scientific << -sten(i+1,j,k,dir::BE) << '\n';

					os << std::setw(width) << (k+1-2)*(nx*ny) + (j+1-2)*nx + i+1-2 << " "
					   << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
					   << std::setw(width) << k + 1 << ", " << " S, "
					   << std::scientific << -sten(i,j,k+1,dir::BNE) << ", "
					   << std::scientific << -sten(i,j,k+1,dir::BN) << ", "
					   << std::scientific << -sten(i+1,j,k+1,dir::BNW) << ", "
					   << std::scientific << -sten(i,j,k,dir::PSW) << ", "
					   << std::scientific << -sten(i,j,k,dir::PS) << ", "
					   << std::scientific << -sten(i+1,j,k,dir::PNW) << ", "
					   << std::scientific << -sten(i,j,k,dir::BSW) << ", "
					   << std::scientific << -sten(i,j,k,dir::BS) << ", "
					   << std::scientific << -sten(i+1,j,k,dir::BSE) << '\n';
				}
			}
		}
	}

	return os;
}

}}
