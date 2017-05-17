#ifndef CEDAR_MULTILEVEL_H
#define CEDAR_MULTILEVEL_H

#include "cedar/types.h"
#include "cedar/util/timer.h"
#include "cedar/cycle/types.h"

namespace cedar {


	template <template<class> class level_container,
		      template<class> class stencil_op,
		      class grid_func,class registry, class fsten, class child>
class multilevel
{
public:
	using conf_ptr = std::shared_ptr<config::reader>;
multilevel(stencil_op<fsten> & fop) : levels(fop), conf(std::make_shared<config::reader>("config.json")) {}
multilevel(stencil_op<fsten> & fop, conf_ptr cfg): levels(fop), conf(cfg) {}
	virtual ~multilevel() {}

	std::shared_ptr<registry> kernel_registry()
	{
		return kreg;
	}

	virtual std::size_t nlevels() { return levels.size(); }

	virtual void log_residual(int lvl, const grid_func & res)
	{
		if (log::info.active()) {
			log::info << "Level " << (levels.size() - lvl - 1) << " residual norm: "
			<< res.template lp_norm<2>() << std::endl;
		}
	}


	virtual void setup_cg_solve()
	{
		auto & cop = levels.get(levels.size() - 1).A;
		kreg->setup_cg_lu(cop, ABD);
		auto kernels = kernel_registry();
		coarse_solver = [&, kernels](grid_func &x, const grid_func & b) {
			kernels->solve_cg(x, b, ABD, bbd);
		};
	}


	virtual void setup_interp(int lvl)
	{
		auto & P = levels.get(lvl+1).P;
		auto & cop = levels.get(lvl+1).A;
		if (lvl == 0) {
			auto & fop = levels.get<fsten>(lvl).A;
			// TODO: check this (lvl changed)
			kreg->setup_interp(lvl+1, lvl, levels.size(), fop, cop, P);
		} else {
			auto & fop = levels.get(lvl).A;
			// TODO: check this (lvl changed)
			kreg->setup_interp(lvl+1, lvl, levels.size(), fop, cop, P);
		}
	}


	virtual void setup_operator(int lvl)
	{
		auto & P = levels.get(lvl+1).P;
		auto & cop = levels.get(lvl+1).A;

		if (lvl == 0) {
			auto & fop = levels.get<fsten>(lvl).A;
			// TODO: check this (lvl changed)
			kreg->galerkin_prod(lvl+1, lvl, levels.size(), P, fop, cop);
		} else {
			auto & fop = levels.get(lvl).A;
			// TODO: check this (lvl changed)
			kreg->galerkin_prod(lvl+1, lvl, levels.size(), P, fop, cop);
		}
	}

	template<class sten>
		void setup_relax_helper(level_container<fsten>::value_type<sten> & level, int lvl)
	{

		auto & sop = level.A;

		std::string relax_type = conf->get<std::string>("solver.relaxation", "point");
		int nrelax_pre = conf->get<int>("solver.cycle.nrelax-pre", 2);
		int nrelax_post = conf->get<int>("solver.cycle.nrelax-post", 1);
		auto kernels = kernel_registry();

		if (relax_type == "point")
			kernels->setup_relax(sop, level.SOR[0]);
		else if (relax_type == "line-x")
			kernels->setup_relax_x(sop, level.SOR[0]);
		else if (relax_type == "line-y")
			kernels->setup_relax_y(sop, level.SOR[0]);
		else if (relax_type == "line-xy") {
			kernels->setup_relax_x(sop, level.SOR[0]);
			kernels->setup_relax_y(sop, level.SOR[1]);
		}
		else if (relax_type == "plane") {
			// TODO: fix this
			//setup_relax_plane(sop, level(lvl));
			// kernels->setup_relax_xy(sop, level(lvl).planes);
		}
		else
			log::error << "Invalid relaxation: " << relax_type << std::endl;

		level.presmoother = [&,lvl,nrelax_pre,kernels,relax_type](const decltype(sop) &A,
		                                                          grid_func &x, const grid_func&b) {
			for (auto i : range(nrelax_pre)) {
				(void) i;
				if (relax_type == "point")
					kernels->relax(A, x, b, level.SOR[0], cycle::Dir::DOWN);
				else if (relax_type == "line-x")
					kernels->relax_lines_x(A, x, b, level.SOR[0], level.res, cycle::Dir::DOWN);
				else if (relax_type == "line-y")
					kernels->relax_lines_y(A, x, b, level.SOR[0], level.res, cycle::Dir::DOWN);
				else if (relax_type == "line-xy") {
					kernels->relax_lines_x(A, x, b, level.SOR[0], level.res, cycle::Dir::DOWN);
					kernels->relax_lines_y(A, x, b, level.SOR[1], level.res, cycle::Dir::DOWN);
				}
				else if (relax_type == "plane") {
					// TODO: fix this
					//relax_plane(A, x, b, cycle::Dir::DOWN, level(lvl));
				}
				else
					log::error << "Invalid relaxation: " << relax_type << std::endl;
			}
		};
		level.postsmoother = [&,lvl,nrelax_post,kernels,relax_type](const decltype(sop) &A,
		                                                            grid_func &x, const grid_func&b) {
			for (auto i: range(nrelax_post)) {
				(void) i;
				if (relax_type == "point")
					kernels->relax(A, x, b, level.SOR[0], cycle::Dir::UP);
				else if (relax_type == "line-x")
					kernels->relax_lines_x(A, x, b, level.SOR[0], level.res, cycle::Dir::UP);
				else if (relax_type == "line-y")
					kernels->relax_lines_y(A, x, b, level.SOR[0], level.res, cycle::Dir::UP);
				else if (relax_type == "line-xy") {
					kernels->relax_lines_y(A, x, b, level.SOR[1], level.res, cycle::Dir::UP);
					kernels->relax_lines_x(A, x, b, level.SOR[0], level.res, cycle::Dir::UP);
				}
				else if (relax_type == "plane") {
					// TODO: fix this
					//relax_plane(av, x, b, cycle::Dir::UP, level(lvl));
				}
				else
					log::error << "Invalid relaxation: " << relax_type << std::endl;
			}
		};
	}

	virtual void setup_relax(int lvl)
	{
		if (lvl == 0) {
			auto & level = levels.get<fsten>(lvl);
			setup_relax_helper(level, lvl);
		} else {
			auto & level = levels.get(lvl);
			setup_relax_helper(level, lvl);
		}
	}


	void setup_space(int nlevels)
	{
		static_cast<child*>(this)->setup_space(nlevels);
	}


	virtual void setup()
	{
		auto num_levels = compute_num_levels(levels.get<fsten>(0).A);
		auto nlevels_conf = conf->get<int>("solver.num-levels", -1);
		if (nlevels_conf > 0) {
			if (nlevels_conf > num_levels) {
				log::error << "too many levels specified" << std::endl;
			} else {
				num_levels = nlevels_conf;
			}
		}
		log::debug << "Using a " << num_levels << " level heirarchy" << std::endl;
		setup_space(num_levels);
		timer_begin("setup");
		for (auto i : range(num_levels-1)) {
			setup_interp(i);
			setup_operator(i);
			setup_relax(i);
		}
		setup_cg_solve();
		timer_end("setup");
	}


	template<class sten>
		void ncycle_helper(level_container<fsten>::value_type<sten> & level,
		                   int lvl, grid_func & x, const grid_func & b, int n)
	{
		auto & A = level.A;

		timer_begin("relaxation");
		level.presmoother(A, x, b);
		timer_end("relaxation");

		grid_func & residual = level.res;
		timer_begin("residual");
		kreg->residual(A, x, b, residual);
		timer_end("residual");

		log_residual(lvl, residual);

		auto & coarse_b = levels.get(lvl+1).b;
		auto & coarse_x = levels.get(lvl+1).x;
		timer_begin("restrict");
		kreg->matvec(levels.get(lvl+1).R, residual, coarse_b);
		timer_end("restrict");
		coarse_x.set(0.0);

		timer_down();

		int coarse_lvl = levels.size() - 1;
		if (lvl+1 == coarse_lvl) {
			timer_begin("coarse-solve");
			coarse_solver(coarse_x, coarse_b);
			timer_end("coarse-solve");
		} else {
			for (auto i : range(n)) {
				(void)i;
				ncycle(lvl+1, coarse_x, coarse_b,n);
			}
		}

		timer_up();

		timer_begin("interp-add");
		//x += levels[lvl].P * coarse_x;
		kreg->interp_add(levels.get(lvl+1).P, coarse_x, level.res, x);
		timer_end("interp-add");

		timer_begin("relaxation");
		level.postsmoother(A, x, b);
		timer_end("relaxation");

		if (log::info.active()) {
			kreg->residual(A, x, b, residual);
			log_residual(lvl, residual);
		}
	}

	void ncycle(int lvl, grid_func & x, const grid_func & b,
		int n=1)
	{
		if (lvl == 0) {
			auto & level = levels.get<fsten>(lvl);
			ncycle_helper(level, lvl, x, b, n);
		} else {
			auto & level = levels.get(lvl);
			ncycle_helper(level, lvl, x, b, n);
		}
	}


	virtual grid_func solve(const grid_func & b)
	{
		auto & level = levels.get<fsten>(0);
		grid_func x = grid_func::zeros_like(b);
		int maxiter = conf->get<int>("solver.max-iter", 10);
		real_t tol = conf->get<real_t>("solver.tol", 1e-8);
		kreg->residual(level.A, x,b,level.res);
		real_t res0_l2 = level.res.template lp_norm<2>();
		log::info << "Initial residual l2 norm: " << res0_l2 << std::endl;

		timer_begin("solve");
		for (auto i: range(maxiter)) {
			vcycle(x, b);
			kreg->residual(level.A, x,b,level.res);
			real_t res_l2 = level.res.template lp_norm<2>();
			real_t rel_l2 = res_l2 / res0_l2;
			log::status << "Iteration " << i << " relative l2 norm: " << rel_l2 << std::endl;
			if (rel_l2 < tol) break;
		}
		timer_end("solve");

		return x;
	}


	virtual void solve(const grid_func & b, grid_func & x)
	{
		auto & level = levels.get<fsten>(0);
		int maxiter = conf->get<int>("solver.max-iter", 10);
		real_t tol = conf->get<real_t>("solver.tol", 1e-8);
		kreg->residual(level.A,x,b,level.res);
		real_t res0_l2 = level.res.template lp_norm<2>();
		log::info << "Initial residual l2 norm: " << res0_l2 << std::endl;

		timer_begin("solve");

		for (auto i: range(maxiter)) {
			vcycle(x, b);
			kreg->residual(level.A,x,b,level.res);
			real_t res_l2 = level.res.template lp_norm<2>();
			real_t rel_l2 = res_l2 / res0_l2;
			log::status << "Iteration " << i << " relative l2 norm: " << rel_l2 << std::endl;
			if (rel_l2 < tol) break;
		}
		timer_end("solve");
	}

	virtual void vcycle(grid_func & x, const grid_func & b)
	{
		if (levels.size() == 1)
			coarse_solver(x, b);
		else
			ncycle(0, x, b);
	}

	int compute_num_levels(stencil_op<fsten> & fop)
	{
		return static_cast<child*>(this)->compute_num_levels(fop);
	}


	/* virtual void setup_relax_plane(stencil_op & sop, LevelType & level) {} */


	/* virtual void relax_plane(const stencil_op & so, */
	/*                          grid_func & x, */
	/*                          const grid_func & b, */
	/*                          cycle::Dir cdir, */
	/*                          LevelType & level) {} */


	config::reader & get_config() { return *conf; }

protected:
	level_container<fsten> levels;
	std::function<void(grid_func &x, const grid_func &b)> coarse_solver;
	std::shared_ptr<config::reader> conf;
	std::shared_ptr<registry> kreg;
	grid_func ABD;
	real_t *bbd;
};

}

#endif
