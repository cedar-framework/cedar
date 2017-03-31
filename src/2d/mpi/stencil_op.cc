#include <cedar/2d/mpi/stencil_op.h>

using namespace cedar;
using namespace cedar::cdr2::mpi;

stencil_op::stencil_op(topo_ptr grd) :
	stencil_op_base(grd->nlocal(0)-2, grd->nlocal(1)-2),
	par_object(grd, grd->comm)
{
	gs.reshape(gs.len(0)+1, gs.len(1)+1, gs.len(2));
	gs.len(0)--;
	gs.len(1)--;
}


grid_func stencil_op::residual(const grid_func &x, const grid_func &b) const
{
	auto r = grid_func(x.grid_ptr());
	stencil_op_base::residual<stencil_op>(x,b,r);

	return r;
}


namespace cedar { namespace cdr2 { namespace mpi {

std::ostream & operator<< (std::ostream & os, const stencil_op &op)
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
				os << std::setw(width) << (jGs+j-2)*NGx + iGs+i-2 << " "
				   << std::setw(width) << iGs + i << ", " << std::setw(width) << jGs + j << ", "
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

}}}
