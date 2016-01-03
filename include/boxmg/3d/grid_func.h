#ifndef BOXMG_3D_GRIDFUNC_H
#define BOXMG_3D_GRIDFUNC_H

#include <iostream>

#include <boxmg/types.h>
#include <boxmg/array.h>
#include <boxmg/grid_quantity.h>

namespace boxmg { namespace bmg3 {

class grid_func : public array<len_t, real_t, 3>, public grid_quantity<len_t, 3>
{
public:
	grid_func(len_t nx, len_t ny, len_t nz, unsigned int nghosts=1);
	grid_func() {}
	static grid_func ones(len_t nx, len_t ny, len_t nz);
	static grid_func zeros(len_t nx, len_t ny, len_t nz);
	static grid_func ones_like(const grid_func & likeable);
	static grid_func zeros_like(const grid_func & likeable);
	friend std::ostream & operator<<(std::ostream &os, const grid_func &obj);
	template<int p> real_t lp_norm() const;
	virtual real_t inf_norm() const;
	grid_func & operator=(grid_func&& gf);
	grid_func(const grid_func & gf) = default;
	grid_func & operator=(const grid_func & gf) = default;
	grid_func & operator -=(const grid_func & rhs);
	friend grid_func operator-(grid_func lhs, const grid_func &rhs) { return lhs -= rhs; }
};


template<int p> real_t grid_func::lp_norm() const
{
	real_t result = 0;

	for (auto &v : vec) result += std::pow(v, p);

	return std::pow(result, 1./p);
}
}}

#endif
