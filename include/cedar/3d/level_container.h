#ifndef CEDAR_3D_LEVEL_CONTAINER_H
#define CEDAR_3D_LEVEL_CONTAINER_H

#include <vector>
#include <cedar/3d/stencil_op.h>

namespace cedar { namespace cdr3 {

	template<template<class> class level, class fsten, class rsten>
	struct get_helper
{
	static level<rsten> & get(std::vector<level<xxvii_pt>> & lvls_xxvii,
	                          std::vector<level<seven_pt>> & lvls_seven,
	                          std::size_t i);
};


template<template<class> class level, class fsten>
	struct get_helper<level, fsten, seven_pt>
{
	static level<seven_pt> & get(std::vector<level<xxvii_pt>> &,
	                            std::vector<level<seven_pt>> & lvls_seven,
	                            std::size_t i)
	{
		if (i != 0) log::error << "coarse operators are nine point (not five)!" << std::endl;
		#ifdef BOUNDS_CHECK
		return lvls_seven.at(0);
		#else
		return lvls_seven[0];
		#endif
	}
};


template<template<class> class level>
	struct get_helper<level, seven_pt, xxvii_pt>
{
	static level<xxvii_pt> & get(std::vector<level<xxvii_pt>> & lvls_xxvii,
	                            std::vector<level<seven_pt>> & lvls_seven,
	                            std::size_t i)
	{
		if (i==0) log::error << "fine grid operator is five point (not nine)!" << std::endl;
		#ifdef BOUNDS_CHECK
		return lvls_xxvii.at(i - lvls_seven.size());
		#else
		return lvls_xxvii[i - lvls_seven.size()];
		#endif
	}
};


template<template<class> class level>
	struct get_helper<level, xxvii_pt, xxvii_pt>
{
	static level<xxvii_pt> & get(std::vector<level<xxvii_pt>> & lvls_xxvii,
	                            std::vector<level<seven_pt>> &,
	                            std::size_t i)
	{
		#ifdef BOUNDS_CHECK
		return lvls_xxvii.at(i);
		#else
		return lvls_xxvii[i];
		#endif
	}
};


template<template<class> class stencil_op, template<class> class level, class sten>
	struct init_helper
	{
		static void init(stencil_op<sten> & fine_op,
		                 std::vector<level<xxvii_pt>> & lvls_xxvii,
		                 std::vector<level<seven_pt>> & lvls_seven,
		                 std::size_t nlevels);
	};


template<template<class> class stencil_op, template<class> class level>
	struct init_helper<stencil_op, level, seven_pt>
	{
		static void init(stencil_op<seven_pt> & fine_op,
		                 std::vector<level<xxvii_pt>> & lvls_xxvii,
		                 std::vector<level<seven_pt>> & lvls_seven,
		                 std::size_t nlevels)
			{
				lvls_seven.emplace_back(fine_op);
				lvls_xxvii.reserve(nlevels-1);
			}
	};


template<template<class> class stencil_op, template<class> class level>
	struct init_helper<stencil_op, level, xxvii_pt>
	{
		static void init(stencil_op<xxvii_pt> & fine_op,
		                 std::vector<level<xxvii_pt>> & lvls_xxvii,
		                 std::vector<level<seven_pt>> &,
		                 std::size_t nlevels)
			{
				lvls_xxvii.reserve(nlevels);
				lvls_xxvii.emplace_back(fine_op);
			}
	};


template<template<class> class level, class fsten>
class level_container
{
public:
	template<class stencil>
		using stencil_op = typename level<fsten>::template stencil_op<stencil>;

	template<class rsten> using level_t = level<rsten>;

level_container(stencil_op<fsten> & fine_op) : fine_op(fine_op) {}

	void init(std::size_t nlevels)
	{
		init_helper<stencil_op, level, fsten>::init(fine_op,
		                                            lvls_xxvii,
		                                            lvls_seven,
		                                            nlevels);
	}

	template<class... Args>
		void add(Args&&... args) {
		lvls_xxvii.emplace_back(std::forward<Args>(args)...);
	}

	template<class rsten=xxvii_pt> level<rsten>&  get(std::size_t i)
	{
		return get_helper<level,fsten, rsten>::get(lvls_xxvii, lvls_seven, i);
	}

	std::size_t size() { return lvls_xxvii.size() + lvls_seven.size(); }

protected:
	stencil_op<fsten> & fine_op;
	std::vector<level<xxvii_pt>> lvls_xxvii;
	std::vector<level<seven_pt>> lvls_seven;
};

}}
#endif
