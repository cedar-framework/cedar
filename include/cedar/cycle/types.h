#ifndef CEDAR_CYCLE_TYPES_H
#define CEDAR_CYCLE_TYPES_H

#include <cedar/types.h>
#include <cedar/kernel_registry.h>
#include <cedar/kernel_manager.h>

namespace cedar { namespace cycle {

template<exec_mode emode, class level_container, class fsten>
class cycle
{
public:
	template<class sten>
		using level_t = typename level_container::template level_t<sten>;
	using grid_func = typename level_t<fsten>::grid_func;
	using stypes = typename level_t<fsten>::stypes;
	using manager = kernel_manager<klist<stypes, emode>, stypes>;

	cycle(level_container & levels,
	      std::function<void(grid_func & x, const grid_func & b)> & coarse_solver) :
		levels(levels),
		coarse_solver(coarse_solver){}
	void set_kernels(std::shared_ptr<manager> kregis) {
		this->kman = kregis;
	}
	virtual void run(grid_func & x, const grid_func & b) = 0;
	void log_residual(std::size_t lvl, const grid_func & res)
	{
		if (log::info.active()) {
			log::info << "Level " << (levels.size() - lvl - 1) << " residual norm: "
			          << const_cast<grid_func&>(res).template lp_norm<2>() << std::endl;
		}
	}

protected:
	std::shared_ptr<manager> kman;
	level_container & levels;
	std::function<void(grid_func & x, const grid_func & b)> & coarse_solver;
};

}}

#endif
