#include <cedar/2d/gpu/prolong_op.h>

using namespace cedar::cdr2::gpu::mpi;

prolong_op::prolong_op(topo_ptr topo) : cdr2::stencil_op<inter_dir>(topo->nlocal(0)-2,
                                                                    topo->nlocal(1)-2),
	par_object(topo, topo->comm)
{
	grid_ = topo;
}


namespace cedar { namespace cdr2 { namespace gpu { namespace mpi {

std::ostream & operator<< (std::ostream &os, const prolong_op &P)
{
	auto & topo = P.grid();
	auto iGs = topo.is(0);
	auto jGs = topo.is(1);
	auto NGx = topo.nglobal(0);

	os << std::setprecision(7);
	for (auto j: P.range(1)) {
		for (auto i: P.range(0)) {
			os << (jGs+j)*NGx + iGs+i << " "
			   << iGs + i << ", " <<  jGs + j << ", "
			   << std::scientific << P(i,j+1,inter_dir::SE) << ", "
			   << std::scientific << P(i,j+1,inter_dir::B)  << ", "
			   << std::scientific << P(i+1,j+1,inter_dir::SW) << ", "
			   << std::scientific << P(i,j,inter_dir::R) << ", "
			   << std::scientific << 1.0 << ", "
			   << std::scientific << P(i+1,j,inter_dir::L) << ", "
			   << std::scientific << P(i,j,inter_dir::NE) << ", "
			   << std::scientific << P(i,j,inter_dir::A) << ", "
			   << std::scientific << P(i+1,j,inter_dir::NW) << '\n';
		}
	}

	return os;
}

}}}}
