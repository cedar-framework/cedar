#ifndef CEDAR_SOLVER_TYPES_H
#define CEDAR_SOLVER_TYPES_H

namespace cedar
{

/**
   Base class for a structured list of types used in a multilevel
   solver.
*/
template<template<class> class stencil_operator,
	class comp_stencil,
	class full_stencil,
	class grid_function,
	class prolong_operator,
	class restrict_operator,
	class relaxation_stencil>
	struct solver_types
	{
		template<class sten>
		using stencil_op = stencil_operator<sten>;
		using comp_sten = comp_stencil;
		using full_sten = full_stencil;
		using grid_func = grid_function;
		using prolong_op = prolong_operator;
		using restrict_op = restrict_operator;
		using relax_stencil = relaxation_stencil;
	};

}

#endif
