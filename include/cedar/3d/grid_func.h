#ifndef CEDAR_3D_GRIDFUNC_H
#define CEDAR_3D_GRIDFUNC_H

#include <iostream>

#include <cedar/types.h>
#include <cedar/array.h>
#include <cedar/grid_quantity.h>

namespace cedar { namespace cdr3 {

class grid_func : public array<real_t, 3>, public grid_quantity<len_t, 3>
{
public:
	grid_func(real_t *ext_data, len_t nx, len_t ny, len_t nz, unsigned int nghosts=1);
	grid_func(len_t nx, len_t ny, len_t nz, unsigned int nghosts=1);
	grid_func() {}
	static grid_func like(const grid_func & likeable);
	static grid_func ones(len_t nx, len_t ny, len_t nz);
	static grid_func zeros(len_t nx, len_t ny, len_t nz);
	static grid_func ones_like(const grid_func & likeable);
	static grid_func zeros_like(const grid_func & likeable);
	friend std::ostream & operator<<(std::ostream &os, const grid_func &obj);
	template<int p> real_t lp_norm() const;
	virtual real_t inf_norm() const;
	grid_func & operator -=(const grid_func & rhs);
	friend grid_func operator-(grid_func lhs, const grid_func &rhs) { return lhs -= rhs; }
};


template<int p> real_t grid_func::lp_norm() const
{
	real_t result = 0;

	for (auto k : this->range(2)) {
		for (auto j : this->range(1)) {
			for (auto i : this->range(0)) {
				result += std::pow((*this)(i,j,k), p);
			}
		}
	}

	return std::pow(result, 1./p);
}
}}

#endif
