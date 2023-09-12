#include <cedar/2d/gpu/stencil_op.h>

namespace cedar { namespace cdr2 { namespace gpu { namespace mpi {

template<>
std::ostream & operator<<(std::ostream & os, const stencil_op<five_pt> & so)
{
	auto & topo = so.grid();
	auto iGs = topo.is(0);
	auto jGs = topo.is(1);
	auto NGx = topo.nglobal(0);
	unsigned int width = 4;

	os << std::setprecision(7);

	for (auto j: so.range(1)) {
		for (auto i: so.range(0)) {
			os << std::setw(width) << (jGs+j-2)*NGx + iGs+i-2 << " "
			   << std::setw(width) << iGs + i << ", " << std::setw(width) << jGs + j << ", "
			   << std::scientific << -so(i  ,j+1,five_pt::s)
			   << std::scientific << -so(i  ,j  ,five_pt::w) << " "
			   << std::scientific <<  so(i  ,j  ,five_pt::c)
			   << std::scientific << -so(i+1,j  ,five_pt::w)
			   << std::scientific << -so(i  ,j  ,five_pt::s) << '\n';
		}
	}

	return os;
}

template<>
std::ostream & operator<<(std::ostream & os, const stencil_op<nine_pt> & so)
{
	auto & topo = so.grid();
	auto iGs = topo.is(0);
	auto jGs = topo.is(1);
	auto NGx = topo.nglobal(0);
	unsigned int width = 4;

	os << std::setprecision(7);

	for (auto j: so.range(1)) {
		for (auto i: so.range(0)) {
			os << std::setw(width) << (jGs+j-2)*NGx + iGs+i-2 << " "
			   << std::setw(width) << iGs + i << ", " << std::setw(width) << jGs + j << ", "
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

}}}}
