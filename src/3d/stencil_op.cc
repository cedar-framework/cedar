#include <iomanip>
#include <cedar/3d/stencil_op.h>


namespace cedar {

template<>
std::ostream & operator<<(std::ostream & os, const cdr3::stencil_op<cdr3::seven_pt> & so)
{
	using namespace cdr3;
	auto nx = so.len(0);
	auto ny = so.len(1);

	unsigned int width = 4;

	os << std::setprecision(7);
	for (auto k : so.range(2)) {
		for (auto j : so.range(1)) {
			for (auto i : so.range(0)) {
				os << std::setw(width) << (k+1-2)*(nx*ny) + (j+1-2)*nx + i+1-2 << " "
				   << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
				   << std::setw(width) << k + 1 << ", "
				   << std::scientific << -so(i,j+1,k,seven_pt::ps) << ", "
				   << std::scientific << -so(i,j,k+1,seven_pt::b) << ", "
				   << std::scientific << -so(i,j,k,seven_pt::pw) << ",  "
				   << std::scientific << so(i,j,k,seven_pt::p) << ", "
				   << std::scientific << -so(i+1,j,k,seven_pt::pw) << ", "
				   << std::scientific << -so(i,j,k,seven_pt::b) << ", "
				   << std::scientific << -so(i,j,k,seven_pt::ps) << '\n';
			}
		}
	}

	return os;
}


template<>
std::ostream & operator<<(std::ostream & os, const cdr3::stencil_op<cdr3::xxvii_pt> & so)
{
	using namespace cdr3;
	auto nx = so.len(0);
	auto ny = so.len(1);

	unsigned int width = 4;

	os << std::setprecision(7);

	for (auto k : so.range(2)) {
		for (auto j : so.range(1)) {
			for (auto i : so.range(0)) {
				os << std::setw(width) << (k+1-2)*(nx*ny) + (j+1-2)*nx + i+1-2 << " "
				   << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
				   << std::setw(width) << k + 1 << ", " << " N, "
				   << std::scientific << -so(i,j+1,k+1,xxvii_pt::bse) << ", "
				   << std::scientific << -so(i,j+1,k+1,xxvii_pt::bs) << ", "
				   << std::scientific << -so(i+1,j+1,k+1,xxvii_pt::bsw) << ", "
				   << std::scientific << -so(i,j+1,k,xxvii_pt::pnw) << ", "
				   << std::scientific << -so(i,j+1,k,xxvii_pt::ps) << ", "
				   << std::scientific << -so(i+1,j+1,k,xxvii_pt::psw) << ", "
				   << std::scientific << -so(i,j+1,k,xxvii_pt::bnw) << ", "
				   << std::scientific << -so(i,j+1,k,xxvii_pt::bn) << ", "
				   << std::scientific << -so(i+1,j+1,k,xxvii_pt::bne) << '\n';

				os << std::setw(width) << (k+1-2)*(nx*ny) + (j+1-2)*nx + i+1-2 << " "
				   << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
				   << std::setw(width) << k + 1 << ", " << " O, "
				   << std::scientific << -so(i,j,k+1,xxvii_pt::be) << ", "
				   << std::scientific << -so(i,j,k+1,xxvii_pt::b) << ", "
				   << std::scientific << -so(i+1,j,k+1,xxvii_pt::bw) << ", "
				   << std::scientific << -so(i,j,k,xxvii_pt::pw) << ", "
				   << std::scientific << so(i,j,k,xxvii_pt::p) << ", "
				   << std::scientific << -so(i+1,j,k,xxvii_pt::pw) << ", "
				   << std::scientific << -so(i,j,k,xxvii_pt::bw) << ", "
				   << std::scientific << -so(i,j,k,xxvii_pt::b) << ", "
				   << std::scientific << -so(i+1,j,k,xxvii_pt::be) << '\n';

				os << std::setw(width) << (k+1-2)*(nx*ny) + (j+1-2)*nx + i+1-2 << " "
				   << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
				   << std::setw(width) << k + 1 << ", " << " S, "
				   << std::scientific << -so(i,j,k+1,xxvii_pt::bne) << ", "
				   << std::scientific << -so(i,j,k+1,xxvii_pt::bn) << ", "
				   << std::scientific << -so(i+1,j,k+1,xxvii_pt::bnw) << ", "
				   << std::scientific << -so(i,j,k,xxvii_pt::psw) << ", "
				   << std::scientific << -so(i,j,k,xxvii_pt::ps) << ", "
				   << std::scientific << -so(i+1,j,k,xxvii_pt::pnw) << ", "
				   << std::scientific << -so(i,j,k,xxvii_pt::bsw) << ", "
				   << std::scientific << -so(i,j,k,xxvii_pt::bs) << ", "
				   << std::scientific << -so(i+1,j,k,xxvii_pt::bse) << '\n';
			}
		}
	}

	return os;
}

}
