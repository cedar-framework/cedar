#ifndef BOXMG_2D_CORE_GRIDFUNC_H
#define BOXMG_2D_CORE_GRIDFUNC_H

#include <cmath>
#include "boxmg-common.h"
#include "array.h"

namespace boxmg { namespace bmg2d { namespace inter {
class prolong_op;
}}}


namespace boxmg { namespace bmg2d {

	class grid_func : public Array<len_t, real_t>
	{
	public:
		using iadd_t = std::tuple<const boxmg::bmg2d::inter::prolong_op&, const grid_func&, const grid_func&>;
		template<typename... ArgTypes>
			grid_func(ArgTypes&&... args) : Array<len_t,real_t>(std::forward<decltype(args)>(args)...){}
		static grid_func ones(len_t nx, len_t ny);
		static grid_func zeros(len_t nx, len_t ny);
		static grid_func zeros_like(const grid_func &likeable);
		static grid_func ones_like(const grid_func &likeable);
		virtual real_t inf_norm() const;
		template<int p> real_t lp_norm() const;
		grid_func & operator+=(iadd_t interp_add_package);
		grid_func & operator-=(const grid_func & rhs);
		friend grid_func operator-(grid_func lhs, const grid_func &rhs) { return lhs -= rhs; }
	};


	template<int p> real_t grid_func::lp_norm() const
	{
		real_t result = 0;

		for (auto &v : data_) result += std::pow(v, p);

		return std::pow(result, 1./p);
	}

}}

#endif
