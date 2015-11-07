#include <iomanip>

#include "inter/types.h"
#include "inter/prolong_op.h"


using namespace boxmg::bmg2d::inter;

ProlongOp::ProlongOp(len_t nx, len_t ny) : StencilOp(nx,ny,true)
{
}

namespace boxmg { namespace bmg2d { namespace inter {
std::ostream & operator<< (std::ostream &os, const ProlongOp &P)
{
	auto & sten = P.stencil();
	// unsigned int width = 8;

	// for (auto j: sten.range(1)) {
	// 	for (auto i: sten.range(0)) {
	// 		os << i << ", " <<  j
	// 		   << std::setw(width) << sten(i,j+1,Dir::SE) << " "
	// 		   << std::setw(width) << sten(i,j+1,Dir::B)  << " "
	// 		   << std::setw(width) << sten(i+1,j+1,Dir::SW) << " "
	// 		   << std::setw(width) << sten(i,j,Dir::R) << " "
	// 		   << std::setw(width) << 1 << " "
	// 		   << std::setw(width) << sten(i+1,j,Dir::L) << " "
	// 		   << std::setw(width) << sten(i,j,Dir::NE) << " "
	// 		   << std::setw(width) << sten(i,j,Dir::A) << " "
	// 		   << std::setw(width) << sten(i+1,j,Dir::NW) << '\n';
	// 	}
	// }

	os << std::setprecision(7);
	for (auto j: sten.range(1)) {
		for (auto i: sten.range(0)) {
			os << i+1 << ", " <<  j+1 << ", "
			   << std::scientific << sten(i,j+1,Dir::SE) << ", "
			   << std::scientific << sten(i,j+1,Dir::B)  << ", "
			   << std::scientific << sten(i+1,j+1,Dir::SW) << ", "
			   << std::scientific << sten(i,j,Dir::R) << ", "
			   << std::scientific << 1.0 << ", "
			   << std::scientific << sten(i+1,j,Dir::L) << ", "
			   << std::scientific << sten(i,j,Dir::NE) << ", "
			   << std::scientific << sten(i,j,Dir::A) << ", "
			   << std::scientific << sten(i+1,j,Dir::NW) << '\n';
		}
	}

	return os;
}


iadd_pack operator*(const ProlongOp & P, const GridFunc & coarse)
{
	return std::make_tuple<std::reference_wrapper<const ProlongOp>,
	                       std::reference_wrapper<const GridFunc>,
	                       std::reference_wrapper<const GridFunc>>(P,coarse,*(P.residual));
}


}}}
