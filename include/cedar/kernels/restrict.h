#ifndef CEDAR_RESTRICT_H
#define CEDAR_RESTRICT_H

#include <cedar/kernel.h>

namespace cedar { namespace kernels {

	template<class solver_types>
	class restriction : public kernel<solver_types>
	{
	public:
		using grid_func = typename kernel<solver_types>::grid_func;
		using restrict_op = typename kernel<solver_types>::restrict_op;

		const std::string name = "restriction";

		virtual void run(const restrict_op & R,
		                 const grid_func & x,
		                 grid_func & y) = 0;
	};
}}

#endif
