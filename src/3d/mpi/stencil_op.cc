#include <cedar/3d/mpi/stencil_op.h>

using namespace cedar;
using namespace cedar::cdr3::mpi;

stencil_op::stencil_op(topo_ptr grd) :
	stencil_op_base(grd->nlocal(0)-2, grd->nlocal(1)-2, grd->nlocal(2)-2),
	par_object(grd, grd->comm)
{
	gs.reshape(gs.len(0)+1, gs.len(1)+1, gs.len(2)+1, gs.len(3));
	gs.len(0)--;
	gs.len(1)--;
	gs.len(2)--;
}


grid_func stencil_op::residual(const grid_func &x, const grid_func &b) const
{
	auto r = grid_func(x.grid_ptr());
	stencil_op_base::residual<stencil_op>(x,b,r);

	return r;
}


namespace cedar { namespace cdr3 { namespace mpi {

std::ostream & operator<< (std::ostream & os, const stencil_op &op)
{
	auto & sten = op.stencil();
	auto & topo = op.grid();
	auto iGs = topo.is(0);
	auto jGs = topo.is(1);
	auto kGs = topo.is(2);
	auto NGx = topo.nglobal(0);
	auto NGy = topo.nglobal(1);
	unsigned int width = 4;

	os << std::setprecision(7);

	if (sten.five_pt()) {
		for (auto k: sten.range(2)) {
			for (auto j: sten.range(1)) {
				for (auto i: sten.range(0)) {
					os << std::setw(width) << (kGs+k-2)*NGx*NGy + (jGs+j-2)*NGx + iGs+i-1 << ", "
					   << std::setw(width) << iGs + i << ", " << std::setw(width) << jGs + j << ", "
					   << std::setw(width) << kGs + k << ", "
					   << std::scientific << -sten(i,j+1,k,dir::PS) << ", "
					   << std::scientific << -sten(i,j,k+1,dir::B) << ", "
					   << std::scientific << -sten(i,j,k,dir::PW) << ", "
					   << std::scientific <<  sten(i,j,k,dir::P) << ", "
					   << std::scientific << -sten(i+1,j,k,dir::PW) << ", "
					   << std::scientific << -sten(i,j,k,dir::B) << ", "
					   << std::scientific << -sten(i,j,k,dir::PS) << '\n';
				}
			}
		}
	} else {
		for (auto k: sten.range(2)) {
			for (auto j: sten.range(1)) {
				for (auto i: sten.range(0)) {
					os << std::setw(width) << (kGs+k-2)*NGx*NGy + (jGs+j-2)*NGx + iGs+i-2 << " "
					   << std::setw(width) << iGs + i << ", " << std::setw(width) << jGs + j << ", "
					   << std::setw(width) << kGs + k << ", N, "
					   << std::scientific << -sten(i,j+1,k+1,dir::BSE) << ", "
					   << std::scientific << -sten(i,j+1,k+1,dir::BS) << ", "
					   << std::scientific << -sten(i+1,j+1,k+1,dir::BSW) << ", "
					   << std::scientific << -sten(i,j+1,k,dir::PNW) << ", "
					   << std::scientific << -sten(i,j+1,k,dir::PS) << ", "
					   << std::scientific << -sten(i+1,j+1,k,dir::PSW) << ", "
					   << std::scientific << -sten(i,j+1,k,dir::BNW) << ", "
					   << std::scientific << -sten(i,j+1,k,dir::BN) << ", "
					   << std::scientific << -sten(i+1,j+1,k,dir::BNE) << '\n';

					os << std::setw(width) << (kGs+k-2)*NGx*NGy + (jGs+j-2)*NGx + iGs+i-2 << " "
					   << std::setw(width) << iGs + i << ", " << std::setw(width) << jGs + j << ", "
					   << std::setw(width) << kGs + k << ", O, "
					   << std::scientific << -sten(i,j,k+1,dir::BE) << ", "
					   << std::scientific << -sten(i,j,k+1,dir::B) << ", "
					   << std::scientific << -sten(i+1,j,k+1,dir::BW) << ", "
					   << std::scientific << -sten(i,j,k,dir::PW) << ", "
					   << std::scientific << sten(i,j,k,dir::P) << ", "
					   << std::scientific << -sten(i+1,j,k,dir::PW) << ", "
					   << std::scientific << -sten(i,j,k,dir::BW) << ", "
					   << std::scientific << -sten(i,j,k,dir::B) << ", "
					   << std::scientific << -sten(i+1,j,k,dir::BE) << '\n';

					os << std::setw(width) << (kGs+k-2)*NGx*NGy + (jGs+j-2)*NGx + iGs+i-2 << " "
					   << std::setw(width) << iGs + i << ", " << std::setw(width) << jGs + j << ", "
					   << std::setw(width) << kGs + k << ", S, "
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

}}}
