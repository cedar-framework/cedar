#include <algorithm>

#include "kernel/name.h"
#include "inter/prolong_op.h"

#include "grid_func.h"

using namespace boxmg::bmg2d::core;


GridFunc GridFunc::ones(len_t nx, len_t ny)
{
	GridFunc ret(nx,ny);

	for (auto j: ret.grange(1)) {
		for (auto i: ret.grange(0)) {
			ret(i,j) = 1.0;
		}
	}

	return ret;
}


GridFunc GridFunc::zeros(len_t nx, len_t ny)
{
	GridFunc ret(nx,ny);

	for (auto j: ret.grange(1)) {
		for (auto i: ret.grange(0)) {
			ret(i,j) = 0.0;
		}
	}

	return ret;
}


GridFunc GridFunc::zeros_like(const GridFunc &like)
{
	GridFunc ret(like.shape(0), like.shape(1));

	for (auto j: ret.grange(1)) {
		for (auto i: ret.grange(0)) {
			ret(i,j) = 0.0;
		}
	}

	return ret;
}


GridFunc GridFunc::ones_like(const GridFunc &like)
{
	GridFunc ret(like.shape(0), like.shape(1));

	for (auto j: ret.grange(1)) {
		for (auto i: ret.grange(0)) {
			ret(i,j) = 1.0;
		}
	}

	return ret;
}


GridFunc & GridFunc::operator+=(iadd_t package)
{
	auto kernels = std::get<0>(package).get_registry();
	kernels->run(kernel::name::interp_add,
	             std::get<0>(package),
	             std::get<1>(package),
	             std::get<2>(package),
	             static_cast<GridFunc&>(*this));

	return *this;
}


boxmg::real_t GridFunc::inf_norm() const
{
	auto abs_compare = [](real_t a, real_t b){
		return (std::abs(a) < std::abs(b));
	};

	auto res = std::max_element(data_.begin(), data_.end(), abs_compare);
	return *res;
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
