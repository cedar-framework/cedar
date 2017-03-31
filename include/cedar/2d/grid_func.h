#ifndef CEDAR_2D_CORE_GRIDFUNC_H
#define CEDAR_2D_CORE_GRIDFUNC_H

#include <cmath>
#include <iostream>

#include <cedar/types.h>
#include <cedar/2d/types.h>
#include <cedar/array.h>
#include <cedar/grid_quantity.h>


namespace cedar { namespace cdr2 { namespace inter {
class prolong_op;
}}}


namespace cedar { namespace cdr2 {

	class grid_func : public array<len_t, real_t, 2>, public grid_quantity<len_t, 2>
	{
	public:
		using iadd_t = std::tuple<const cedar::cdr2::inter::prolong_op&, const grid_func&, const grid_func&>;
		using array<len_t, real_t, 2>::operator();
		grid_func(len_t nx, len_t ny, unsigned int nghosts=1);
		grid_func() {}
		grid_func & operator=(grid_func&& gf);
		grid_func(const grid_func & gf) = default;
		grid_func & operator=(const grid_func & gf) = default;
		static grid_func ones(len_t nx, len_t ny);
		static grid_func zeros(len_t nx, len_t ny);
		static grid_func random(len_t nx, len_t ny);
		static grid_func like(const grid_func &likeable);
		static grid_func zeros_like(const grid_func &likeable);
		static grid_func ones_like(const grid_func &likeable);
		virtual real_t inf_norm() const;
		template<int p> real_t lp_norm() const;
		grid_func & operator+=(iadd_t interp_add_package);
		grid_func & operator-=(const grid_func & rhs);
		friend grid_func operator-(grid_func lhs, const grid_func &rhs) { return lhs -= rhs; }
		friend std::ostream & operator<<(std::ostream &os, const grid_func & obj);
	};


	template<int p> real_t grid_func::lp_norm() const
	{
		real_t result = 0;

		for (auto j : this->range(1)) {
			for (auto i : this->range(0)) {
				result += std::pow((*this)(i,j), p);
			}
		}

		return std::pow(result, 1./p);
	}

}}

#endif