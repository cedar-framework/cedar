#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/inter/mpi/prolong_op.h>

using namespace cedar::cdr2::mpi;

grid_func::grid_func(topo_ptr grid) :
	::cedar::cdr2::grid_func(grid->nlocal(0)-2,grid->nlocal(1)-2), par_object(grid, grid->comm) {}

grid_func::grid_func(len_t nx, len_t ny): ::cedar::cdr2::grid_func(nx,ny)
{}


grid_func grid_func::like(const grid_func &likeable)
{
	grid_func ret(likeable.grid_ptr());

	return ret;
}


grid_func grid_func::zeros_like(const grid_func &like)
{
	grid_func ret(like.grid_ptr());

	for (auto j: ret.grange(1)) {
		for (auto i: ret.grange(0)) {
			ret(i,j) = 0.0;
		}
	}

	return ret;
}


grid_func grid_func::ones_like(const grid_func &like)
{
	grid_func ret(like.grid_ptr());

	for (auto j: ret.grange(1)) {
		for (auto i: ret.grange(0)) {
			ret(i,j) = 1.0;
		}
	}

	return ret;
}


grid_func & grid_func::operator-=(const grid_func &rhs)
{
	for (auto j: this->range(1)) {
		for (auto i: this->range(0)) {
			(*this)(i,j) -= rhs(i,j);
		}
	}

	return *this;
}


cedar::real_t grid_func::inf_norm() const
{
	auto mval = cdr2::grid_func::inf_norm();

	MPI_Allreduce(MPI_IN_PLACE, &mval, 1, MPI_DOUBLE, MPI_MAX, grid_->comm);

	return mval;
}

grid_func & grid_func::operator+=(iadd_t package)
{
	auto kernels = std::get<0>(package).get_registry();
	kernels->run(kernel_name::interp_add,
	             std::get<0>(package),
	             std::get<1>(package),
	             std::get<2>(package),
	             static_cast<grid_func&>(*this));

	return *this;
}

namespace cedar { namespace cdr2 { namespace mpi {
std::ostream & operator<< (std::ostream &os, const grid_func & obj)
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
