#ifndef CEDAR_CYCLE_FCYCLE_H
#define CEDAR_CYCLE_FCYCLE_H

#include <cedar/types.h>
#include <cedar/util/timer.h>
#include <cedar/cycle/vcycle.h>

namespace cedar { namespace cycle {

template<exec_mode emode, class level_container, class fsten>
class fcycle : public cycle<emode, level_container, fsten>
{
	template<class sten>
		using level_t = typename level_container::template level_t<sten>;
	template<class sten>
		using stencil_op = typename level_t<fsten>::template stencil_op<sten>;
	using grid_func = typename level_t<fsten>::grid_func;
	// extract kernel names for convenience
	using stypes = typename level_t<fsten>::stypes;

	using coarsen_op = kernels::coarsen_op<stypes>;
	using interp_add = kernels::interp_add<stypes>;
	using residual = kernels::residual<stypes>;
	using restriction = kernels::restriction<stypes>;
	using setup_interp = kernels::setup_interp<stypes>;
	using solve_cg = kernels::solve_cg<stypes>;

	using parent = cycle<emode, level_container, fsten>;
	using parent::levels;
	using parent::coarse_solver;
	using parent::log_residual;
	using parent::kman;

public:
fcycle(level_container & levels, std::function<void(grid_func&,const grid_func&)> & coarse_solver) :
parent::cycle(levels, coarse_solver), vcycler(levels, coarse_solver) {
}
	virtual void run(grid_func & x, const grid_func & b) override
	{
		vcycler.set_kernels(kman);

		if (levels.size() == 1)
			coarse_solver(x, b);
		else
			fmg_cycle(0, x, b);
	}


	void fmg_cycle(std::size_t lvl, grid_func & x, const grid_func & b)
	{
		if (lvl == levels.size() - 1) {
			// /\* x.set(0.0); *\/
			coarse_solver(x, b);
			if (log::info.active()) {
				auto & A = levels.get(levels.size() - 1).A;
				auto & res = levels.get(levels.size()-1).res;
				kman->template run<residual>(A, x, b, res);
				log_residual(lvl, res);
			}
		} else {
			if (lvl == 0) {
				auto & level = levels.template get<fsten>(lvl);
				fmg_cycle_helper(lvl, level, x, b);
			} else {
				auto & level = levels.get(lvl);
				fmg_cycle_helper(lvl, level, x, b);
			}
		}
	}


	template<class sten>
		void fmg_cycle_helper(std::size_t lvl, level_t<sten> & level, grid_func & x, const grid_func & b)
	{
			auto & coarse_b = levels.get(lvl+1).b;
			auto & coarse_x = levels.get(lvl+1).x;
			kman->template run<restriction>(levels.get(lvl+1).R, b, coarse_b);
			fmg_cycle(lvl+1, coarse_x, coarse_b);
			x.set(0.0);
			level.res.set(0.0);
			kman->template run<interp_add>(levels.get(lvl+1).P, coarse_x, level.res, x);
			vcycler.ncycle(lvl, x, b);
	}

private:
	vcycle<emode, level_container, fsten> vcycler;

};

}}

#endif
