#include "grid_func.h"

using namespace boxmg::bmg2d::mpi;

GridFunc::GridFunc(topo_ptr grid) :
	::boxmg::bmg2d::GridFunc(grid->nlocal(0)-2,grid->nlocal(1)-2), comm(grid->comm), grid_(grid) {}


GridFunc GridFunc::zeros_like(const GridFunc &like)
{
	GridFunc ret(like.grid_ptr());

	for (auto j: ret.grange(1)) {
		for (auto i: ret.grange(0)) {
			ret(i,j) = 0.0;
		}
	}

	return ret;
}


GridFunc GridFunc::ones_like(const GridFunc &like)
{
	GridFunc ret(like.grid_ptr());

	for (auto j: ret.grange(1)) {
		for (auto i: ret.grange(0)) {
			ret(i,j) = 1.0;
		}
	}

	return ret;
}


GridFunc & GridFunc::operator-=(const GridFunc &rhs)
{
	for (auto j: this->range(1)) {
		for (auto i: this->range(0)) {
			(*this)(i,j) -= rhs(i,j);
		}
	}

	return *this;
}


boxmg::real_t GridFunc::inf_norm() const
{
	auto mval = bmg2d::GridFunc::inf_norm();

	MPI_Allreduce(MPI_IN_PLACE, &mval, 1, MPI_DOUBLE, MPI_MAX, grid_->comm);

	return mval;
}

namespace boxmg { namespace bmg2d { namespace mpi {
std::ostream & operator<< (std::ostream &os, const GridFunc & obj)
{
	auto & topo = obj.grid();
	auto iGs = topo.is(0);
	auto jGs = topo.is(1);
	auto NGx = topo.nglobal(0);
	unsigned int width = 4;

	os << std::setprecision(7);

	for (auto j: obj.range(1)) {
		for (auto i: obj.range(0)) {
			os << std::setw(width) << (jGs+j-2)*NGx + iGs+i-2 << " "
			   << std::setw(width) << iGs + i << ", " << std::setw(width) << jGs + j << ", "
			   << std::scientific << obj(i,j) << '\n';
		}
	}

	return os;
}
}}}
