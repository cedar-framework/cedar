#include <boxmg/2d/inter/mpi/prolong_op.h>

using namespace boxmg::bmg2d::inter::mpi;

prolong_op::prolong_op(bmg2d::mpi::topo_ptr topo) : inter::prolong_op(topo->nlocal(0)-2,
                                                             topo->nlocal(1)-2),
                                            grid_(topo)
{
}




namespace boxmg { namespace bmg2d { namespace inter { namespace mpi {

std::ostream & operator<< (std::ostream &os, const prolong_op &P)
{
	auto & topo = P.grid();
	auto iGs = topo.is(0);
	auto jGs = topo.is(1);
	auto NGx = topo.nglobal(0);

	auto & sten = P.stencil();
	// unsigned int width = 8;

	// for (auto j: sten.range(1)) {
	// 	for (auto i: sten.range(0)) {
	// 		os << i << ", " <<  j
	// 		   << std::setw(width) << sten(i,j+1,Dir::SE) << " "
	// 		   << std::setw(width) << sten(i,j+1,Dir::B)  << " "
	// 		   << std::setw(width) << sten(i+1,j+1,Dir::SW) << " "
	// 		   << std::setw(width) << sten(i,j,Dir::R) << " "
	// 		   << std::setw(width) << 1 << " "
	// 		   << std::setw(width) << sten(i+1,j,Dir::L) << " "
	// 		   << std::setw(width) << sten(i,j,Dir::NE) << " "
	// 		   << std::setw(width) << sten(i,j,Dir::A) << " "
	// 		   << std::setw(width) << sten(i+1,j,Dir::NW) << '\n';
	// 	}
	// }

	os << std::setprecision(7);
	for (auto j: sten.range(1)) {
		for (auto i: sten.range(0)) {
			os << (jGs+j)*NGx + iGs+i << " "
			   << iGs + i << ", " <<  jGs + j << ", "
			   << std::scientific << sten(i,j+1,dir::SE) << ", "
			   << std::scientific << sten(i,j+1,dir::B)  << ", "
			   << std::scientific << sten(i+1,j+1,dir::SW) << ", "
			   << std::scientific << sten(i,j,dir::R) << ", "
			   << std::scientific << 1.0 << ", "
			   << std::scientific << sten(i+1,j,dir::L) << ", "
			   << std::scientific << sten(i,j,dir::NE) << ", "
			   << std::scientific << sten(i,j,dir::A) << ", "
			   << std::scientific << sten(i+1,j,dir::NW) << '\n';
		}
	}

	return os;
}

}}}}
