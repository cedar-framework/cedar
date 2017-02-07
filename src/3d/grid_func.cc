#include <boxmg/kernel_name.h>
#include <boxmg/3d/inter/prolong_op.h>

#include <boxmg/3d/grid_func.h>


using namespace boxmg;
using namespace boxmg::bmg3;


grid_func & grid_func::operator=(grid_func &&gf)
{
	vec = std::move(gf.vec);
	strides = std::move(gf.strides);
	extents = std::move(gf.extents);
	num_ghosts = gf.num_ghosts;
	range_ = std::move(gf.range_);
	grange_ = std::move(gf.grange_);

	return *this;
}


grid_func::grid_func(len_t nx, len_t ny, len_t nz, unsigned int nghosts) :
	array<len_t, real_t, 3>(nx+2*nghosts, ny+2*nghosts, nz+2*nghosts)
{
	num_ghosts = nghosts;
	range_[0] = ::boxmg::range(static_cast<len_t>(nghosts), static_cast<len_t>(nx + nghosts));
	range_[1] = ::boxmg::range(static_cast<len_t>(nghosts), static_cast<len_t>(ny + nghosts));
	range_[2] = ::boxmg::range(static_cast<len_t>(nghosts), static_cast<len_t>(nz + nghosts));

	grange_[0] = ::boxmg::range(static_cast<len_t>(0), static_cast<len_t>(nx + 3*nghosts));
	grange_[1] = ::boxmg::range(static_cast<len_t>(0), static_cast<len_t>(ny + 3*nghosts));
	grange_[2] = ::boxmg::range(static_cast<len_t>(0), static_cast<len_t>(nz + 3*nghosts));
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


boxmg::real_t grid_func::inf_norm() const
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


namespace boxmg { namespace bmg3 {
std::ostream & operator<<(std::ostream &os, const grid_func &obj)
{
	for (auto k : obj.range(2)) {
		for (auto j : obj.range(1)) {
			for (auto i : obj.range(0)) {
				os << i << " " << j << " " << k << " => " << std::to_string(obj(i,j,k)) << '\n';
			}
		}
	}

	return os;
}
}}
