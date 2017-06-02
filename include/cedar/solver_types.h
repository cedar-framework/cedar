#ifndef CEDAR_SOLVER_TYPES_H
#define CEDAR_SOLVER_TYPES_H

namespace cedar
{

template<template<class> class stencil_operator,
	class grid_function,
	class prolong_operator,
	class restrict_operator,
	class relaxation_stencil>
	struct solver_types
	{
		template<class sten>
		using stencil_op = stencil_operator<sten>;
		using grid_func = grid_function;
		using prolong_op = prolong_operator;
		using restrict_op = restrict_operator;
		using relax_stencil = relaxation_stencil;
	};

}

#endif
