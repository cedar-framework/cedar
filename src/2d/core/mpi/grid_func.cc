#include "grid_func.h"

using namespace boxmg::bmg2d;
using namespace boxmg::bmg2d::core::mpi;

GridFunc::GridFunc(topo_ptr grid) :
	::boxmg::bmg2d::core::GridFunc(grid->nlocal(0)-2,grid->nlocal(1)-2), comm(grid->comm), grid_(grid) {}


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
	auto mval = core::GridFunc::inf_norm();

	MPI_Allreduce(MPI_IN_PLACE, &mval, 1, MPI_DOUBLE, MPI_MAX, grid_->comm);

	return mval;
}
