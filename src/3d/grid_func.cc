#include <cedar/3d/grid_func.h>


using namespace cedar;
using namespace cedar::cdr3;


grid_func::grid_func(len_t nx, len_t ny, len_t nz, unsigned int nghosts) :
	array<real_t, 3>(nx+2*nghosts, ny+2*nghosts, nz+2*nghosts)
{
	num_ghosts = nghosts;
	range_[0] = ::cedar::range(static_cast<len_t>(nghosts), static_cast<len_t>(nx + nghosts));
	range_[1] = ::cedar::range(static_cast<len_t>(nghosts), static_cast<len_t>(ny + nghosts));
	range_[2] = ::cedar::range(static_cast<len_t>(nghosts), static_cast<len_t>(nz + nghosts));

	grange_[0] = ::cedar::range(static_cast<len_t>(0), static_cast<len_t>(nx + 2*nghosts));
	grange_[1] = ::cedar::range(static_cast<len_t>(0), static_cast<len_t>(ny + 2*nghosts));
	grange_[2] = ::cedar::range(static_cast<len_t>(0), static_cast<len_t>(nz + 2*nghosts));
}


grid_func::grid_func(real_t *ext_data, len_t nx, len_t ny, len_t nz, unsigned int nghosts) :
	array<real_t, 3>(ext_data, nx+2*nghosts, ny+2*nghosts, nz+2*nghosts)
{
	num_ghosts = nghosts;
	range_[0] = ::cedar::range(static_cast<len_t>(nghosts), static_cast<len_t>(nx + nghosts));
	range_[1] = ::cedar::range(static_cast<len_t>(nghosts), static_cast<len_t>(ny + nghosts));
	range_[2] = ::cedar::range(static_cast<len_t>(nghosts), static_cast<len_t>(nz + nghosts));

	grange_[0] = ::cedar::range(static_cast<len_t>(0), static_cast<len_t>(nx + 2*nghosts));
	grange_[1] = ::cedar::range(static_cast<len_t>(0), static_cast<len_t>(ny + 2*nghosts));
	grange_[2] = ::cedar::range(static_cast<len_t>(0), static_cast<len_t>(nz + 2*nghosts));
}


grid_func grid_func::ones(len_t nx, len_t ny, len_t nz)
{
	grid_func ret(nx, ny, nz);

	ret.set(1.0);

	return ret;
}


grid_func grid_func::like(const grid_func &likeable)
{
	grid_func ret(likeable.shape(0), likeable.shape(1), likeable.shape(2));

	return ret;
}


grid_func grid_func::zeros(len_t nx, len_t ny, len_t nz)
{
	grid_func ret(nx, ny, nz);

	ret.set(0.0);

	return ret;
}


grid_func grid_func::zeros_like(const grid_func &like)
{
	grid_func ret(like.shape(0), like.shape(1), like.shape(2));

	ret.set(0.0);

	return ret;
}


grid_func grid_func::ones_like(const grid_func &like)
{
	grid_func ret(like.shape(0), like.shape(1), like.shape(2));

	ret.set(1.0);

	return ret;
}


cedar::real_t grid_func::inf_norm() const
{
	auto abs_compare = [](real_t a, real_t b){
		return (std::abs(a) < std::abs(b));
	};

	real_t cmax = 0;
	for (auto k : this->range(2)) {
		for (auto j : this->range(1)) {
			for (auto i : this->range(0)) {
				if (abs_compare(cmax, (*this)(i,j,k)))
					cmax = (*this)(i,j,k);
			}
		}
	}

	return cmax;
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

namespace cedar { namespace cdr3 {
std::ostream & operator<<(std::ostream &os, const grid_func &obj)
{
	unsigned int width = 4;

	os << std::setprecision(7);

	for (auto k : obj.range(2)) {
		for (auto j : obj.range(1)) {
			for (auto i : obj.range(0)) {
				// os << i << " " << j << " " << k << " => " << std::to_string(obj(i,j,k)) << '\n';
				os << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
				   << std::setw(width) << k+1 << ", "
				   << std::scientific << obj(i,j,k) << '\n';
			}
		}
	}

	return os;
}
}}
