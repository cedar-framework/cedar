#include <boxmg/3d/inter/mpi/prolong_op.h>

using namespace boxmg::bmg3::inter::mpi;

prolong_op::prolong_op(topo_ptr topo) : mpi::stencil_op(topo->nlocal(0)-2,
                                                        topo->nlocal(1)-2,
                                                        topo->nlocal(2)-2, true)
{
	grid_ = topo;
}




namespace boxmg { namespace bmg3 { namespace inter { namespace mpi {

std::ostream & operator<< (std::ostream &os, const prolong_op &P)
{
	auto & topo = P.grid();
	auto iGs = topo.is(0);
	auto jGs = topo.is(1);
	auto kGs = topo.is(2);
	auto NGx = topo.nglobal(0);
	auto NGy = topo.nglobal(1);

	auto & sten = P.stencil();

	os << std::setprecision(7);
	for (auto k: sten.range(2)) {
		for (auto j: sten.range(1)) {
			for (auto i: sten.range(0)) {
				os << (kGs+k)*NGx*NGy + (jGs+j)*NGx + iGs+i << " "
				   << iGs + i << ", " <<  jGs + j << ", " << kGs + k << " N, "
				   << std::scientific << sten(i,j,k+1,dir::BNE) << ", "
				   << std::scientific << sten(i,j,k+1,dir::XZSE) << ", "
				   << std::scientific << sten(i,j+1,k+1,dir::BSE) << ", "
				   << std::scientific << sten(i,j,k,dir::XYNE) << ", "
				   << std::scientific << sten(i,j,k,dir::XYR) << ", "
				   << std::scientific << sten(i,j+1,k,dir::XYSE) << ", "
				   << std::scientific << sten(i,j,k,dir::TNE) << ", "
				   << std::scientific << sten(i,j,k,dir::XZNE) << ", "
				   << std::scientific << sten(i,j+1,k,dir::TSE) << '\n';

				os << (kGs+k)*NGx*NGy + (jGs+j)*NGx + iGs+i << " "
				   << iGs + i << ", " <<  jGs + j << ", " << kGs + k << " O, "
				   << std::scientific << sten(i,j,k+1,dir::YZSW) << ", "
				   << std::scientific << sten(i,j,k+1,dir::XZB) << ", "
				   << std::scientific << sten(i,j+1,k+1,dir::YZSE) << ", "
				   << std::scientific << sten(i,j,k,dir::XYA) << ", "
				   << std::scientific << 1.0 << ", "
				   << std::scientific << sten(i,j+1,k,dir::XYB) << ", "
				   << std::scientific << sten(i,j,k,dir::YZNW) << ", "
				   << std::scientific << sten(i,j,k,dir::XZA) << ", "
				   << std::scientific << sten(i,j+1,k,dir::YZNE) << '\n';

				os << (kGs+k)*NGx*NGy + (jGs+j)*NGx + iGs+i << " "
				   << iGs + i << ", " <<  jGs + j << ", " << kGs + k << " S, "
				   << std::scientific << sten(i+1,j,k+1,dir::BNW) << ", "
				   << std::scientific << sten(i+1,j,k+1,dir::XZSW) << ", "
				   << std::scientific << sten(i+1,j+1,k+1,dir::BSW) << ", "
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

iadd_pack operator*(const prolong_op & P, const mpi::grid_func & coarse)
{
	return std::make_tuple<std::reference_wrapper<const prolong_op>,
	                       std::reference_wrapper<const mpi::grid_func>,
	                       std::reference_wrapper<const mpi::grid_func>>(P,coarse,*(P.residual));
}

}}}}
