#include <cedar/3d/mpi/prolong_op.h>

using namespace cedar::cdr3::mpi;

prolong_op::prolong_op(topo_ptr topo) : cdr3::stencil_op<inter_dir>(topo->nlocal(0)-2,
                                                                    topo->nlocal(1)-2,
                                                                    topo->nlocal(2)-2),
	par_object(topo, topo->comm)
{
	grid_ = topo;
}




namespace cedar { namespace cdr3 { namespace mpi {

std::ostream & operator<< (std::ostream &os, const prolong_op &P)
{
	auto & topo = P.grid();
	auto iGs = topo.is(0);
	auto jGs = topo.is(1);
	auto kGs = topo.is(2);
	auto NGx = topo.nglobal(0);
	auto NGy = topo.nglobal(1);


	os << std::setprecision(7);
	for (auto k: P.range(2)) {
		for (auto j: P.range(1)) {
			for (auto i: P.range(0)) {
				os << (kGs+k)*NGx*NGy + (jGs+j)*NGx + iGs+i << " "
				   << iGs + i << ", " <<  jGs + j << ", " << kGs + k << " N, "
				   << std::scientific << P(i,j,k+1,inter_dir::BNE) << ", "
				   << std::scientific << P(i,j,k+1,inter_dir::XZSE) << ", "
				   << std::scientific << P(i,j+1,k+1,inter_dir::BSE) << ", "
				   << std::scientific << P(i,j,k,inter_dir::XYNE) << ", "
				   << std::scientific << P(i,j,k,inter_dir::XYR) << ", "
				   << std::scientific << P(i,j+1,k,inter_dir::XYSE) << ", "
				   << std::scientific << P(i,j,k,inter_dir::TNE) << ", "
				   << std::scientific << P(i,j,k,inter_dir::XZNE) << ", "
				   << std::scientific << P(i,j+1,k,inter_dir::TSE) << '\n';

				os << (kGs+k)*NGx*NGy + (jGs+j)*NGx + iGs+i << " "
				   << iGs + i << ", " <<  jGs + j << ", " << kGs + k << " O, "
				   << std::scientific << P(i,j,k+1,inter_dir::YZSW) << ", "
				   << std::scientific << P(i,j,k+1,inter_dir::XZB) << ", "
				   << std::scientific << P(i,j+1,k+1,inter_dir::YZSE) << ", "
				   << std::scientific << P(i,j,k,inter_dir::XYA) << ", "
				   << std::scientific << 1.0 << ", "
				   << std::scientific << P(i,j+1,k,inter_dir::XYB) << ", "
				   << std::scientific << P(i,j,k,inter_dir::YZNW) << ", "
				   << std::scientific << P(i,j,k,inter_dir::XZA) << ", "
				   << std::scientific << P(i,j+1,k,inter_dir::YZNE) << '\n';

				os << (kGs+k)*NGx*NGy + (jGs+j)*NGx + iGs+i << " "
				   << iGs + i << ", " <<  jGs + j << ", " << kGs + k << " S, "
				   << std::scientific << P(i+1,j,k+1,inter_dir::BNW) << ", "
				   << std::scientific << P(i+1,j,k+1,inter_dir::XZSW) << ", "
				   << std::scientific << P(i+1,j+1,k+1,inter_dir::BSW) << ", "
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

}}}
