#include <cedar/3d/mpi/stencil_op.h>

using namespace cedar;
using namespace cedar::cdr3::mpi;

namespace cedar { namespace cdr3 { namespace mpi {

template<>
std::ostream & operator<< (std::ostream & os, const stencil_op<seven_pt> &op)
{
	auto & topo = op.grid();
	auto iGs = topo.is(0);
	auto jGs = topo.is(1);
	auto kGs = topo.is(2);
	auto NGx = topo.nglobal(0);
	auto NGy = topo.nglobal(1);
	unsigned int width = 4;

	os << std::setprecision(7);

	for (auto k: op.range(2)) {
		for (auto j: op.range(1)) {
			for (auto i: op.range(0)) {
				os << std::setw(width) << (kGs+k-2)*NGx*NGy + (jGs+j-2)*NGx + iGs+i-1 << ", "
				   << std::setw(width) << iGs + i << ", " << std::setw(width) << jGs + j << ", "
				   << std::setw(width) << kGs + k << ", "
				   << std::scientific << -op(i,j+1,k,seven_pt::ps) << ", "
				   << std::scientific << -op(i,j,k+1,seven_pt::b) << ", "
				   << std::scientific << -op(i,j,k,seven_pt::pw) << ", "
				   << std::scientific <<  op(i,j,k,seven_pt::p) << ", "
				   << std::scientific << -op(i+1,j,k,seven_pt::pw) << ", "
				   << std::scientific << -op(i,j,k,seven_pt::b) << ", "
				   << std::scientific << -op(i,j,k,seven_pt::ps) << '\n';
			}
		}
	}
	return os;
}

template<>
std::ostream & operator<< (std::ostream & os, const stencil_op<xxvii_pt> &op)
{
	auto & topo = op.grid();
	auto iGs = topo.is(0);
	auto jGs = topo.is(1);
	auto kGs = topo.is(2);
	auto NGx = topo.nglobal(0);
	auto NGy = topo.nglobal(1);
	unsigned int width = 4;

	os << std::setprecision(7);

	for (auto k: op.range(2)) {
		for (auto j: op.range(1)) {
			for (auto i: op.range(0)) {
				os << std::setw(width) << (kGs+k-2)*NGx*NGy + (jGs+j-2)*NGx + iGs+i-2 << " "
				   << std::setw(width) << iGs + i << ", " << std::setw(width) << jGs + j << ", "
				   << std::setw(width) << kGs + k << ", N, "
				   << std::scientific << -op(i,j+1,k+1,xxvii_pt::bse) << ", "
				   << std::scientific << -op(i,j+1,k+1,xxvii_pt::bs) << ", "
				   << std::scientific << -op(i+1,j+1,k+1,xxvii_pt::bsw) << ", "
				   << std::scientific << -op(i,j+1,k,xxvii_pt::pnw) << ", "
				   << std::scientific << -op(i,j+1,k,xxvii_pt::ps) << ", "
				   << std::scientific << -op(i+1,j+1,k,xxvii_pt::psw) << ", "
				   << std::scientific << -op(i,j+1,k,xxvii_pt::bnw) << ", "
				   << std::scientific << -op(i,j+1,k,xxvii_pt::bn) << ", "
				   << std::scientific << -op(i+1,j+1,k,xxvii_pt::bne) << '\n';

				os << std::setw(width) << (kGs+k-2)*NGx*NGy + (jGs+j-2)*NGx + iGs+i-2 << " "
				   << std::setw(width) << iGs + i << ", " << std::setw(width) << jGs + j << ", "
				   << std::setw(width) << kGs + k << ", O, "
				   << std::scientific << -op(i,j,k+1,xxvii_pt::be) << ", "
				   << std::scientific << -op(i,j,k+1,xxvii_pt::b) << ", "
				   << std::scientific << -op(i+1,j,k+1,xxvii_pt::bw) << ", "
				   << std::scientific << -op(i,j,k,xxvii_pt::pw) << ", "
				   << std::scientific << op(i,j,k,xxvii_pt::p) << ", "
				   << std::scientific << -op(i+1,j,k,xxvii_pt::pw) << ", "
				   << std::scientific << -op(i,j,k,xxvii_pt::bw) << ", "
				   << std::scientific << -op(i,j,k,xxvii_pt::b) << ", "
				   << std::scientific << -op(i+1,j,k,xxvii_pt::be) << '\n';

				os << std::setw(width) << (kGs+k-2)*NGx*NGy + (jGs+j-2)*NGx + iGs+i-2 << " "
				   << std::setw(width) << iGs + i << ", " << std::setw(width) << jGs + j << ", "
				   << std::setw(width) << kGs + k << ", S, "
				   << std::scientific << -op(i,j,k+1,xxvii_pt::bne) << ", "
				   << std::scientific << -op(i,j,k+1,xxvii_pt::bn) << ", "
				   << std::scientific << -op(i+1,j,k+1,xxvii_pt::bnw) << ", "
				   << std::scientific << -op(i,j,k,xxvii_pt::psw) << ", "
				   << std::scientific << -op(i,j,k,xxvii_pt::ps) << ", "
				   << std::scientific << -op(i+1,j,k,xxvii_pt::pnw) << ", "
				   << std::scientific << -op(i,j,k,xxvii_pt::bsw) << ", "
				   << std::scientific << -op(i,j,k,xxvii_pt::bs) << ", "
				   << std::scientific << -op(i+1,j,k,xxvii_pt::bse) << '\n';
			}
		}
	}

	return os;
}

}}}
