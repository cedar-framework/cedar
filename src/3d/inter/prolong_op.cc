#include <iomanip>

#include <cedar/3d/inter/types.h>
#include <cedar/3d/inter/prolong_op.h>

using namespace cedar::cdr3::inter;


prolong_op::prolong_op(len_t nx, len_t ny, len_t nz) : stencil_op<inter::dir>(nx,ny,nz) {}

namespace cedar { namespace cdr3 { namespace inter {

std::ostream & operator<<(std::ostream & os, const prolong_op & P)
{
	os << std::setprecision(7);
	for (auto k : P.range(2)) {
		for (auto j : P.range(1)) {
			for (auto i : P.range(0)) {
				os << i+1 << ", " << j+1 << ", " << k+1 << " N, "
				   << std::scientific << P(i,j,k+1, dir::BNE) << ", "
				   << std::scientific << P(i,j,k+1, dir::XZSE) << ", "
				   << std::scientific << P(i,j+1,k+1,dir::BSE) << ", "
				   << std::scientific << P(i,j,k,dir::XYNE) << ", "
				   << std::scientific << P(i,j,k,dir::XYR) << ", "
				   << std::scientific << P(i,j+1,k,dir::XYSE) << ", "
				   << std::scientific << P(i,j,k,dir::TNE) << ", "
				   << std::scientific << P(i,j,k,dir::XZNE) << ", "
				   << std::scientific << P(i,j+1,k,dir::TSE) << '\n';

				os << i+1 << ", " << j+1 << ", " << k+1 << " O, "
				   << std::scientific << P(i,j,k+1,dir::YZSW) << ", "
				   << std::scientific << P(i,j,k+1,dir::XZB) << ", "
				   << std::scientific << P(i,j+1,k+1,dir::YZSE) << ", "
				   << std::scientific << P(i,j,k,dir::XYA) << ", "
				   << std::scientific << 1.0 << ", "
				   << std::scientific << P(i,j+1,k,dir::XYB) << ", "
				   << std::scientific << P(i,j,k,dir::YZNW) << ", "
				   << std::scientific << P(i,j,k,dir::XZA) << ", "
				   << std::scientific << P(i,j+1,k,dir::YZNE) << '\n';

				os << i+1 << ", " << j+1 << ", " << k+1 << " S, "
				   << std::scientific << P(i+1,j,k+1,dir::BNW) << ", "
				   << std::scientific << P(i+1,j,k+1,dir::XZSW) << ", "
				   << std::scientific << P(i+1,j+1,k+1, dir::BSW) << ", "
				   << std::scientific << P(i+1,j,k,dir::XYNW) << ", "
				   << std::scientific << P(i+1,j,k,dir::XYL) << ", "
				   << std::scientific << P(i+1,j+1,k,dir::XYSW) << ", "
				   << std::scientific << P(i+1,j,k,dir::TNW) << ", "
				   << std::scientific << P(i+1,j,k,dir::XZNW) << ", "
				   << std::scientific << P(i+1,j+1,k,dir::TSW) << '\n';
			}
		}
	}

	return os;
}

}}}
