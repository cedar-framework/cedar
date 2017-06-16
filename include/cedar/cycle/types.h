#ifndef CEDAR_CYCLE_TYPES_H
#define CEDAR_CYCLE_TYPES_H

#include <cedar/types.h>


namespace cedar { namespace cycle {

enum class Dir {UP, DOWN};

template<class level_container,
	class registry, class fsten>
class cycle
{
	template<class sten>
		using level_t = typename level_container::template level_t<sten>;
	template<class sten>
		using stencil_op = typename level_t<fsten>::template stencil_op<sten>;
	using grid_func = typename level_t<fsten>::grid_func;

public:
cycle(level_container & levels,
      std::function<void(grid_func & x, const grid_func & b)> & coarse_solver) :
		levels(levels),
		coarse_solver(coarse_solver){}
		void set_registry(std::shared_ptr<registry> kregis) {
			this->kreg = kregis;
		}
virtual void run(grid_func & x, const grid_func & b) = 0;
void log_residual(std::size_t lvl, const grid_func & res)
{
	if (log::info.active()) {
		log::info << "Level " << (levels.size() - lvl - 1) << " residual norm: "
		          << res.template lp_norm<2>() << std::endl;
	}
}

protected:
std::shared_ptr<registry> kreg;
	level_container & levels;
	std::function<void(grid_func & x, const grid_func & b)> & coarse_solver;
};

}}

#endif
