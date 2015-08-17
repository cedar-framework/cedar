#ifndef BOXMG_2D_CORE_GRIDFUNC_H
#define BOXMG_2D_CORE_GRIDFUNC_H

#include <cmath>
#include "boxmg-common.h"
#include "array.h"

namespace boxmg { namespace bmg2d { namespace inter {
class ProlongOp;
}}}


namespace boxmg { namespace bmg2d { namespace core {

	class GridFunc : public Array<len_t, real_t>
	{
	public:
		using iadd_t = std::tuple<const boxmg::bmg2d::inter::ProlongOp&, const GridFunc&, const GridFunc&>;
		template<typename... ArgTypes>
			GridFunc(ArgTypes&&... args) : Array<len_t,real_t>(std::forward<decltype(args)>(args)...){}
		static GridFunc ones(len_t nx, len_t ny);
		static GridFunc zeros(len_t nx, len_t ny);
		static GridFunc zeros_like(const GridFunc &likeable);
		static GridFunc ones_like(const GridFunc &likeable);
		virtual real_t inf_norm() const;
		template<int p> real_t lp_norm() const;
		GridFunc & operator+=(iadd_t interp_add_package);
		GridFunc & operator-=(const GridFunc & rhs);
		friend GridFunc operator-(GridFunc lhs, const GridFunc &rhs) { return lhs -= rhs; }
	};


	template<int p> real_t GridFunc::lp_norm() const
	{
		real_t result = 0;

		for (auto &v : data_) result += std::pow(v, p);

		return std::pow(result, 1./p);
	}

}}}

#endif
