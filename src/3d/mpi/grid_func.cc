#include <cedar/3d/mpi/grid_func.h>
#include <cedar/3d/inter/mpi/prolong_op.h>

using namespace cedar::cdr3::mpi;

grid_func::grid_func(topo_ptr grid) :
	::cedar::cdr3::grid_func(grid->nlocal(0)-2, grid->nlocal(1)-2, grid->nlocal(2) - 2),
	par_object(grid, grid->comm) {}


grid_func::grid_func(len_t nx, len_t ny, len_t nz) : ::cedar::cdr3::grid_func(nx,ny,nz)
{}


grid_func grid_func::zeros(topo_ptr grid)
{
	grid_func ret(grid);

	ret.set(0.0);

	return ret;
}


grid_func grid_func::ones(topo_ptr grid)
{
	grid_func ret(grid);

	ret.set(1.0);

	return ret;
}


grid_func grid_func::like(const grid_func &likeable)
{
	grid_func ret(likeable.grid_ptr());

	return ret;
}


grid_func grid_func::zeros_like(const grid_func &like)
{
	grid_func ret(like.grid_ptr());

	ret.set(0.0);

	return ret;
}


grid_func grid_func::ones_like(const grid_func &like)
{
	grid_func ret(like.grid_ptr());

	ret.set(1.0);

	return ret;
}


grid_func & grid_func::operator-=(const grid_func &rhs)
{
	for (auto k : this->range(2)) {
		for (auto j: this->range(1)) {
			for (auto i: this->range(0)) {
				(*this)(i,j,k) -= rhs(i,j,k);
			}
		}
	}

	return *this;
}


cedar::real_t grid_func::inf_norm() const
{
	auto mval = cdr3::grid_func::inf_norm();

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

namespace cedar { namespace cdr3 { namespace mpi {
std::ostream & operator<< (std::ostream &os, const grid_func & obj)
{
	auto & topo = obj.grid();
	auto iGs = topo.is(0);
	auto jGs = topo.is(1);
	auto kGs = topo.is(2);
	auto NGx = topo.nglobal(0);
	auto NGy = topo.nglobal(1);
	unsigned int width = 4;

	os << std::setprecision(7);

	for (auto k : obj.range(2)) {
		for (auto j: obj.range(1)) {
			for (auto i: obj.range(0)) {
				os << std::setw(width) << (kGs+k-2)*NGx*NGy + (jGs+j-2)*NGx + iGs+i-2 << ", "
				   << std::setw(width) << iGs + i << ", " << std::setw(width) << jGs + j << ", "
				   << std::setw(width) << kGs + k << ", "
				   << std::scientific << obj(i,j,k) << '\n';
			}
		}
	}

	return os;
}

		}}}
