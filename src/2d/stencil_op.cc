#include <cedar/2d/stencil_op.h>

namespace cedar
{
	template<>
	std::ostream & operator<<(std::ostream & os, const cdr2::stencil_op<cdr2::five_pt> & so)
	{
		using namespace cdr2;
		unsigned int width = 4;

		os << std::setprecision(7);

		for (auto j: so.range(1)) {
			for (auto i: so.range(0)) {
				os << std::setw(width) << (j+1-2)*so.len(0) + i+1-2 << " "
				   << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
				   << std::scientific << -so(i  ,j+1,five_pt::s)
				   << std::scientific << -so(i  ,j,  five_pt::w) << " "
				   << std::scientific <<  so(i  ,j,  five_pt::c)
				   << std::scientific << -so(i+1,j,  five_pt::w)
				   << std::scientific << -so(i  ,j,  five_pt::s) << '\n';
			}
		}
		return os;
	}

	template<>
	std::ostream & operator<<(std::ostream & os, const cdr2::stencil_op<cdr2::nine_pt> & so)
	{
		using namespace cdr2;
		unsigned int width = 4;

		os << std::setprecision(7);

		for (auto j: so.range(1)) {
			for (auto i: so.range(0)) {
				os << std::setw(width) << (j+1-2)*so.len(0) + i+1-2 << " "
				   << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
				   << std::scientific << -so(i  ,j+1,nine_pt::nw)
				   << std::scientific << -so(i  ,j+1,nine_pt::s)
				   << std::scientific << -so(i+1,j+1,nine_pt::sw)
				   << std::scientific << -so(i  ,j  ,nine_pt::w) << " "
				   << std::scientific <<  so(i  ,j  ,nine_pt::c)
				   << std::scientific << -so(i+1,j  ,nine_pt::w)
				   << std::scientific << -so(i  ,j  ,nine_pt::sw)
				   << std::scientific << -so(i  ,j  ,nine_pt::s)
				   << std::scientific << -so(i+1,j  ,nine_pt::nw) << '\n';
			}
		}
		return os;
	}
}
