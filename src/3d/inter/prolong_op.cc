#include <iomanip>

#include <boxmg/3d/inter/types.h>
#include <boxmg/3d/inter/prolong_op.h>

using namespace boxmg::bmg3::inter;


prolong_op::prolong_op(len_t nx, len_t ny, len_t nz) : stencil_op(nx,ny,nz,true) {}

namespace boxmg { namespace bmg3 { namespace inter {

std::ostream & operator<<(std::ostream & os, const prolong_op & P)
{
	auto & sten = P.stencil();

	os << std::setprecision(7);
	for (auto k : sten.range(2)) {
		for (auto j : sten.range(1)) {
			for (auto i : sten.range(0)) {
				os << i+1 << ", " << j+1 << ", " << k+1 << " N, "
				   << std::scientific << sten(i,j,k+1, dir::BNE) << ", "
				   << std::scientific << sten(i,j,k+1, dir::XZSE) << ", "
				   << std::scientific << sten(i,j+1,k+1,dir::BSE) << ", "
				   << std::scientific << sten(i,j,k,dir::XYNE) << ", "
				   << std::scientific << sten(i,j,k,dir::XYR) << ", "
				   << std::scientific << sten(i,j+1,k,dir::XYSE) << ", "
				   << std::scientific << sten(i,j,k,dir::TNE) << ", "
				   << std::scientific << sten(i,j,k,dir::XZNE) << ", "
				   << std::scientific << sten(i,j+1,k,dir::TSE) << '\n';

				os << i+1 << ", " << j+1 << ", " << k+1 << " O, "
				   << std::scientific << sten(i,j,k+1,dir::YZSW) << ", "
				   << std::scientific << sten(i,j,k+1,dir::XZB) << ", "
				   << std::scientific << sten(i,j+1,k+1,dir::YZSE) << ", "
				   << std::scientific << sten(i,j,k,dir::XYA) << ", "
				   << std::scientific << 1.0 << ", "
				   << std::scientific << sten(i,j+1,k,dir::XYB) << ", "
				   << std::scientific << sten(i,j,k,dir::YZNW) << ", "
				   << std::scientific << sten(i,j,k,dir::XZA) << ", "
				   << std::scientific << sten(i,j+1,k,dir::YZNE) << '\n';

				os << i+1 << ", " << j+1 << ", " << k+1 << " S, "
				   << std::scientific << sten(i+1,j,k+1,dir::BNW) << ", "
				   << std::scientific << sten(i+1,j,k+1,dir::XZSW) << ", "
				   << std::scientific << sten(i+1,j+1,k+1, dir::BSW) << ", "
				   << std::scientific << sten(i+1,j,k,dir::XYNW) << ", "
				   << std::scientific << sten(i+1,j,k,dir::XYL) << ", "
				   << std::scientific << sten(i+1,j+1,k,dir::XYSW) << ", "
				   << std::scientific << sten(i+1,j,k,dir::TNW) << ", "
				   << std::scientific << sten(i+1,j,k,dir::XZNW) << ", "
				   << std::scientific << sten(i+1,j+1,k,dir::TSW) << '\n';
			}
		}
	}

	return os;
}


iadd_pack operator*(const prolong_op & P, const grid_func & coarse)
{
	return std::make_tuple<std::reference_wrapper<const prolong_op>,
	                       std::reference_wrapper<const grid_func>,
	                       std::reference_wrapper<const grid_func>>(P,coarse,*(P.residual));
}

}}}
