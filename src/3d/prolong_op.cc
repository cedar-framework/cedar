#include <iomanip>

#include <cedar/3d/prolong_op.h>

using namespace cedar::cdr3;


prolong_op::prolong_op(len_t nx, len_t ny, len_t nz) : stencil_op<inter_dir>(nx,ny,nz) {}

namespace cedar { namespace cdr3 {

std::ostream & operator<<(std::ostream & os, const prolong_op & P)
{
	os << std::setprecision(7);
	for (auto k : P.range(2)) {
		for (auto j : P.range(1)) {
			for (auto i : P.range(0)) {
				os << i+1 << ", " << j+1 << ", " << k+1 << " N, "
				   << std::scientific << P(i,j,k+1, inter_dir::BNE) << ", "
				   << std::scientific << P(i,j,k+1, inter_dir::XZSE) << ", "
				   << std::scientific << P(i,j+1,k+1,inter_dir::BSE) << ", "
				   << std::scientific << P(i,j,k,inter_dir::XYNE) << ", "
				   << std::scientific << P(i,j,k,inter_dir::XYR) << ", "
				   << std::scientific << P(i,j+1,k,inter_dir::XYSE) << ", "
				   << std::scientific << P(i,j,k,inter_dir::TNE) << ", "
				   << std::scientific << P(i,j,k,inter_dir::XZNE) << ", "
				   << std::scientific << P(i,j+1,k,inter_dir::TSE) << '\n';

				os << i+1 << ", " << j+1 << ", " << k+1 << " O, "
				   << std::scientific << P(i,j,k+1,inter_dir::YZSW) << ", "
				   << std::scientific << P(i,j,k+1,inter_dir::XZB) << ", "
				   << std::scientific << P(i,j+1,k+1,inter_dir::YZSE) << ", "
				   << std::scientific << P(i,j,k,inter_dir::XYA) << ", "
				   << std::scientific << 1.0 << ", "
				   << std::scientific << P(i,j+1,k,inter_dir::XYB) << ", "
				   << std::scientific << P(i,j,k,inter_dir::YZNW) << ", "
				   << std::scientific << P(i,j,k,inter_dir::XZA) << ", "
				   << std::scientific << P(i,j+1,k,inter_dir::YZNE) << '\n';

				os << i+1 << ", " << j+1 << ", " << k+1 << " S, "
				   << std::scientific << P(i+1,j,k+1,inter_dir::BNW) << ", "
				   << std::scientific << P(i+1,j,k+1,inter_dir::XZSW) << ", "
				   << std::scientific << P(i+1,j+1,k+1, inter_dir::BSW) << ", "
				   << std::scientific << P(i+1,j,k,inter_dir::XYNW) << ", "
				   << std::scientific << P(i+1,j,k,inter_dir::XYL) << ", "
				   << std::scientific << P(i+1,j+1,k,inter_dir::XYSW) << ", "
				   << std::scientific << P(i+1,j,k,inter_dir::TNW) << ", "
				   << std::scientific << P(i+1,j,k,inter_dir::XZNW) << ", "
				   << std::scientific << P(i+1,j+1,k,inter_dir::TSW) << '\n';
			}
		}
	}

	return os;
}

}}
