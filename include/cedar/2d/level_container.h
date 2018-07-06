#ifndef CEDAR_2D_LEVEL_CONTAINER_H
#define CEDAR_2D_LEVEL_CONTAINER_H

#include <vector>
#include <cedar/2d/stencil_op.h>
#include <cedar/services/mempool.h>

namespace cedar { namespace cdr2 {

	template<template<class> class level, class fsten, class rsten>
	struct get_helper
{
	static level<rsten> & get(std::vector<level<nine_pt>> & lvls_nine,
	                          std::vector<level<five_pt>> & lvls_five,
	                          std::size_t i);
};


template<template<class> class level, class fsten>
	struct get_helper<level, fsten, five_pt>
{
	static level<five_pt> & get(std::vector<level<nine_pt>> & lvls_nine,
	                            std::vector<level<five_pt>> & lvls_five,
	                            std::size_t i)
	{
		if (i != 0) log::error << "coarse operators are nine point (not five)!" << std::endl;
		#ifdef BOUNDS_CHECK
		return lvls_five.at(0);
		#else
		return lvls_five[0];
		#endif
	}
};


template<template<class> class level>
	struct get_helper<level, five_pt, nine_pt>
{
	static level<nine_pt> & get(std::vector<level<nine_pt>> & lvls_nine,
	                            std::vector<level<five_pt>> & lvls_five,
	                            std::size_t i)
	{
		if (i==0) log::error << "fine grid operator is five point (not nine)!" << std::endl;
		#ifdef BOUNDS_CHECK
		return lvls_nine.at(i - lvls_five.size());
		#else
		return lvls_nine[i - lvls_five.size()];
		#endif
	}
};


template<template<class> class level>
	struct get_helper<level, nine_pt, nine_pt>
{
	static level<nine_pt> & get(std::vector<level<nine_pt>> & lvls_nine,
	                            std::vector<level<five_pt>> & lvls_five,
	                            std::size_t i)
	{
		#ifdef BOUNDS_CHECK
		return lvls_nine.at(i);
		#else
		return lvls_nine[i];
		#endif
	}
};


template<template<class> class stencil_op, template<class> class level, class sten>
	struct init_helper
	{
		static void init(stencil_op<sten> & fine_op,
		                 std::vector<level<nine_pt>> & lvls_nine,
		                 std::vector<level<five_pt>> & lvls_five,
		                 std::size_t nlevels);
	};


template<template<class> class stencil_op, template<class> class level>
	struct init_helper<stencil_op, level, five_pt>
	{
		static void init(stencil_op<five_pt> & fine_op,
		                 services::mempool & mpool,
		                 std::vector<level<nine_pt>> & lvls_nine,
		                 std::vector<level<five_pt>> & lvls_five,
		                 std::size_t nlevels)
			{
				lvls_five.emplace_back(mpool, fine_op);
				lvls_nine.reserve(nlevels-1);
			}
	};


template<template<class> class stencil_op, template<class> class level>
	struct init_helper<stencil_op, level, nine_pt>
	{
		static void init(stencil_op<nine_pt> & fine_op,
		                 services::mempool & mpool,
		                 std::vector<level<nine_pt>> & lvls_nine,
		                 std::vector<level<five_pt>> & lvls_five,
		                 std::size_t nlevels)
			{
				lvls_nine.reserve(nlevels);
				lvls_nine.emplace_back(mpool, fine_op);
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

	void init(services::mempool & mpool, std::size_t nlevels)
	{
		init_helper<stencil_op, level, fsten>::init(fine_op,
		                                            mpool,
		                                            lvls_nine,
		                                            lvls_five,
		                                            nlevels);
	}

	template<class... Args>
		void add(Args&&... args) {
		lvls_nine.emplace_back(std::forward<Args>(args)...);
	}

	template<class rsten=nine_pt> level<rsten>&  get(std::size_t i)
	{
		return get_helper<level,fsten, rsten>::get(lvls_nine, lvls_five, i);
	}

	std::size_t size() { return lvls_nine.size() + lvls_five.size(); }

protected:
	stencil_op<fsten> & fine_op;
	std::vector<level<nine_pt>> lvls_nine;
	std::vector<level<five_pt>> lvls_five;
};

}}
#endif
