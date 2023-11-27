#ifndef CEDAR_CYCLE_VCYCLE_H
#define CEDAR_CYCLE_VCYCLE_H

#include <cedar/types.h>
#include <cedar/util/timer.h>

namespace cedar { namespace cycle {

template<exec_mode emode, class level_container, class fsten>
class vcycle : public cycle<emode, level_container, fsten>
{
public:
	template<class sten>
		using level_t = typename level_container::template level_t<sten>;
	template<class sten>
		using stencil_op = typename level_t<fsten>::template stencil_op<sten>;
	using grid_func = typename level_t<fsten>::grid_func;
	// extract kernel names for convenience
	using stypes = typename level_t<fsten>::stypes;

	using interp_add = kernels::interp_add<stypes>;
	using restriction = kernels::restriction<stypes>;
	using residual = kernels::residual<stypes>;
	using setup_interp = kernels::setup_interp<stypes>;
	using solve_cg = kernels::solve_cg<stypes>;

	using parent = cycle<emode, level_container, fsten>;
	using parent::levels;
	using parent::coarse_solver;
	using parent::log_residual;
	using parent::kman;

	vcycle(level_container & levels, std::function<void(grid_func&,const grid_func&)> & coarse_solver) :
		parent::cycle(levels, coarse_solver) {}
	virtual void run(grid_func & x, const grid_func & b) override
	{
		if (levels.size() == 1)
			coarse_solver(x, b);
		else
			ncycle(0, x, b);
	}


	void ncycle(std::size_t lvl, grid_func & x, const grid_func & b,
		int n=1)
	{
		if (lvl == 0) {
			auto & level = levels.template get<fsten>(lvl);
			ncycle_helper(level, lvl, x, b, n);
		} else {
			auto & level = levels.get(lvl);
			ncycle_helper(level, lvl, x, b, n);
		}
	}


	template<class sten>
		void ncycle_helper(level_t<sten> & level,
		                   std::size_t lvl, grid_func & x, const grid_func & b, int n)
	{
		auto & A = level.A;

                // std::cerr << "Before smoothing: " << std::endl << x.to_buffer() << std::endl;

		timer_begin("relaxation");
		level.presmoother(A, x, b);
		timer_end("relaxation");

		grid_func & res = level.res;
		timer_begin("residual");
		kman->template run<residual>(A, x, b, res);
		timer_end("residual");

                // std::cerr << "Post smoothing: " << std::endl << x.to_buffer() << std::endl;
                // std::cerr << "Residual: " << std::endl << res.to_buffer() << std::endl;

		log_residual(lvl, res);

		auto & coarse_b = levels.get(lvl+1).b;
		auto & coarse_x = levels.get(lvl+1).x;
		timer_begin("restrict");

		kman->template run<restriction>(levels.get(lvl+1).R, res, coarse_b);
		timer_end("restrict");
		coarse_x.set(0.0);

                // std::cerr << "Coarse solution: " << std::endl << coarse_x.to_buffer() << std::endl;
                // std::cerr << "coarse rhs: " << std::endl << coarse_b.to_buffer() << std::endl;

		timer_down();

		std::size_t coarse_lvl = levels.size() - 1;
		if (lvl+1 == coarse_lvl) {
			timer_begin("coarse-solve");
			coarse_solver(coarse_x, coarse_b);
			timer_end("coarse-solve");
			if (log::info.active()) {
				auto & coarse_A = levels.get(levels.size() - 1).A;
				auto & coarse_residual = levels.get(levels.size()-1).res;
				kman->template run<residual>(coarse_A, coarse_x, coarse_b, coarse_residual);
				log_residual(lvl+1, coarse_residual);
			}
		} else {
			for (auto i : range(n)) {
				(void)i;
				ncycle(lvl+1, coarse_x, coarse_b,n);
			}
		}

		timer_up();

                // std::cerr << "Coarse solution after solve: " << std::endl << coarse_x.to_buffer() << std::endl;

		timer_begin("interp-add");
		kman->template run<interp_add>(levels.get(lvl+1).P, coarse_x, res, x);
                // std::cerr << "Interpolated solution: " << std::endl << x.to_buffer() << std::endl;
		timer_end("interp-add");

		timer_begin("relaxation");
		level.postsmoother(A, x, b);
		timer_end("relaxation");

		if (log::info.active()) {
			kman->template run<residual>(A, x, b, res);
			log_residual(lvl, res);
                        std::cerr << "residual: " << std::endl << res.to_buffer() << std::endl;
		}
	}
};


}}


#endif
