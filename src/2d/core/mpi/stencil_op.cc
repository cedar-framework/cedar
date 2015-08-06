#include "stencil_op.h"

using namespace boxmg::bmg2d;
using namespace boxmg::bmg2d::core::mpi;

StencilOp::StencilOp(topo_ptr grd) :
	core::StencilOp(grd->nlocal(0)-2, grd->nlocal(1)-2),
	comm(grd->comm), grid_(grd)
{
	gs.stride(1)++;
	gs.stride(2) = (gs.len(0)+1)*(gs.len(1)+1);
}


GridFunc StencilOp::residual(const GridFunc &x, const GridFunc &b) const
{
	auto r = GridFunc(x.grid_ptr());
	core::StencilOp::residual(x,b,r);

	return r;
}


void StencilOp::residual(const GridFunc &x, const GridFunc &b, GridFunc &r) const
{
	core::StencilOp::residual(x,b,r);
}


namespace boxmg { namespace bmg2d { namespace core { namespace mpi {

std::ostream & operator<< (std::ostream & os, const StencilOp &op)
{
	auto & sten = op.stencil();
	auto & topo = op.grid();
	auto iGs = topo.is(0);
	auto jGs = topo.is(1);
	auto NGx = topo.nglobal(0);
	unsigned int width = 4;

	os << std::setprecision(7);

	if (sten.five_pt()) {
		for (auto j: sten.range(1)) {
			for (auto i: sten.range(0)) {
				os << std::setw(width) << (jGs+j-2)*NGx + iGs+i-2 << " "
				   << std::setw(width) << iGs + i << ", " << std::setw(width) << jGs + j << ", "
				   << std::scientific << -sten(i,j,Dir::N)
				   << std::scientific << -sten(i,j,Dir::W) << " "
				   << std::scientific <<  sten(i,j,Dir::C)
				   << std::scientific << -sten(i,j,Dir::E)
				   << std::scientific << -sten(i,j,Dir::S) << '\n';
			}
		}
	} else {
		for (auto j: sten.range(1)) {
			for (auto i: sten.range(0)) {
				os << std::setw(width) << (jGs+j-2)*NGx + iGs+i-2 << " "
				   << std::setw(width) << iGs + i << ", " << std::setw(width) << jGs + j << ", "
				   << std::scientific << -sten(i,j,Dir::NW)
				   << std::scientific << -sten(i,j,Dir::N)
				   << std::scientific << -sten(i,j,Dir::NE)
				   << std::scientific << -sten(i,j,Dir::W) << " "
				   << std::scientific <<  sten(i,j,Dir::C)
				   << std::scientific << -sten(i,j,Dir::E)
				   << std::scientific << -sten(i,j,Dir::SW)
				   << std::scientific << -sten(i,j,Dir::S)
				   << std::scientific << -sten(i,j,Dir::SE) << '\n';
			}
		}
	}

	return os;
}

}}}}
