#ifndef BOXMG_MULTILEVEL_H
#define BOXMG_MULTILEVEL_H

#include "boxmg/types.h"
#include "boxmg/util/timer.h"
#include "boxmg/cycle/types.h"
#include "boxmg/discrete_op.h"

namespace boxmg {


template <class LevelType,class stencil_op, class grid_func,class registry>
class multilevel
{
public:
	using conf_ptr = std::shared_ptr<config::reader>;
multilevel() : conf(std::make_shared<config::reader>("config.json")) {}
multilevel(conf_ptr cfg): conf(cfg) {}
	virtual ~multilevel() {}

	std::shared_ptr<registry> kernel_registry()
	{
		return kreg;
	}

	virtual LevelType & level(int i)
	{
		int ind;
		if (i < 0) ind = (i + 1)*-1;
		else ind = levels.size() - i - 1;
		#ifdef DEBUG
		return levels.at(ind);
		#else
		return levels[ind];
		#endif
	}
	virtual const LevelType & level(int i) const
	{
		int ind;
		if (i < 0) ind = (i + 1)*-1;
		else ind = levels.size() - i - 1;
		#ifdef DEBUG
		return levels.at(ind);
		#else
		return levels[ind];
		#endif
	}
	virtual int nlevels() { return levels.size(); }

	virtual void log_residual(int lvl, const grid_func & res)
	{
		if (log::info.active()) {
			log::info << "Level " << (levels.size() - lvl - 1) << " residual norm: "
			<< res.template lp_norm<2>() << std::endl;
		}
	}


	virtual void setup_cg_solve()
	{
		kreg->setup_cg_lu(levels.back().A, ABD);
		auto kernels = kernel_registry();
		coarse_solver = [&, kernels](const discrete_op<grid_func> &A, grid_func &x, const grid_func & b) {
			kernels->solve_cg(x, b, ABD, bbd);
			const stencil_op &av = dynamic_cast<const stencil_op&>(A);
			grid_func & residual = levels[levels.size()-1].res;
			av.residual(x,b,residual);
			log::info << "Level 0 residual norm: " << residual.template lp_norm<2>() << std::endl;
		};
	}


	virtual void setup_interp(int lvl)
	{
		auto & P = level(lvl).P;
		auto & fop = level(lvl).A;
		auto & cop = level(lvl-1).A;
		kreg->setup_interp(lvl+1, lvl, levels.size(), fop, cop, P);
	}


	virtual void setup_operator(int lvl)
	{
		auto & P = level(lvl).P;
		auto & fop = level(lvl).A;
		auto & cop = level(lvl-1).A;
		kreg->galerkin_prod(lvl+1, lvl, levels.size(), P, fop, cop);
	}


	virtual void setup_relax(int lvl)
	{
		auto & sop = level(lvl).A;

		std::string relax_type = conf->get<std::string>("solver.relaxation", "point");
		int nrelax_pre = conf->get<int>("solver.cycle.nrelax-pre", 2);
		int nrelax_post = conf->get<int>("solver.cycle.nrelax-post", 1);
		auto kernels = kernel_registry();

		if (relax_type == "point")
			kernels->setup_relax(sop, level(lvl).SOR[0]);
		else if (relax_type == "line-x")
			kernels->setup_relax_x(sop, level(lvl).SOR[0]);
		else if (relax_type == "line-y")
			kernels->setup_relax_y(sop, level(lvl).SOR[0]);
		else if (relax_type == "line-xy") {
			kernels->setup_relax_x(sop, level(lvl).SOR[0]);
			kernels->setup_relax_y(sop, level(lvl).SOR[1]);
		}
		else if (relax_type == "plane") {
			setup_relax_plane(sop, level(lvl));
			// kernels->setup_relax_xy(sop, level(lvl).planes);
		}
		else
			log::error << "Invalid relaxation: " << relax_type << std::endl;

		level(lvl).presmoother = [&,lvl,nrelax_pre,kernels,relax_type](const discrete_op<grid_func> &A, grid_func &x, const grid_func&b) {
			const stencil_op & av = dynamic_cast<const stencil_op &>(A);
			for (auto i : range(nrelax_pre)) {
				(void) i;
				if (relax_type == "point")
					kernels->relax(av, x, b, level(lvl).SOR[0], cycle::Dir::DOWN);
				else if (relax_type == "line-x")
					kernels->relax_lines_x(av, x, b, level(lvl).SOR[0], level(lvl).res, cycle::Dir::DOWN);
				else if (relax_type == "line-y")
					kernels->relax_lines_y(av, x, b, level(lvl).SOR[0], level(lvl).res, cycle::Dir::DOWN);
				else if (relax_type == "line-xy") {
					kernels->relax_lines_x(av, x, b, level(lvl).SOR[0], level(lvl).res, cycle::Dir::DOWN);
					kernels->relax_lines_y(av, x, b, level(lvl).SOR[1], level(lvl).res, cycle::Dir::DOWN);
				}
				else if (relax_type == "plane") {
					relax_plane(av, x, b, cycle::Dir::DOWN, level(lvl));
				}
				else
					log::error << "Invalid relaxation: " << relax_type << std::endl;
			}
		};
		level(lvl).postsmoother = [&,lvl,nrelax_post,kernels,relax_type](const discrete_op<grid_func> &A, grid_func &x, const grid_func&b) {

			const stencil_op & av = dynamic_cast<const stencil_op &>(A);
			for (auto i: range(nrelax_post)) {
				(void) i;
				if (relax_type == "point")
					kernels->relax(av, x, b, level(lvl).SOR[0], cycle::Dir::UP);
				else if (relax_type == "line-x")
					kernels->relax_lines_x(av, x, b, level(lvl).SOR[0], level(lvl).res, cycle::Dir::UP);
				else if (relax_type == "line-y")
					kernels->relax_lines_y(av, x, b, level(lvl).SOR[0], level(lvl).res, cycle::Dir::UP);
				else if (relax_type == "line-xy") {
					kernels->relax_lines_y(av, x, b, level(lvl).SOR[1], level(lvl).res, cycle::Dir::UP);
					kernels->relax_lines_x(av, x, b, level(lvl).SOR[0], level(lvl).res, cycle::Dir::UP);
				}
				else if (relax_type == "plane") {
					relax_plane(av, x, b, cycle::Dir::UP, level(lvl));
				}
				else
					log::error << "Invalid relaxation: " << relax_type << std::endl;
			}
		};
	}


	// eventually make this necessary
	virtual void setup_space(int nlevels)
	{
	}


	virtual void setup(stencil_op&& fop)
	{
		levels.emplace_back(std::move(fop));
		levels.back().A.set_registry(kreg);
		auto num_levels = compute_num_levels(levels[0].A);
		auto nlevels_conf = conf->get<int>("solver.num-levels", -1);
		if (nlevels_conf > 0) {
			if (nlevels_conf > num_levels) {
				log::error << "too many levels specified" << std::endl;
			} else {
				num_levels = nlevels_conf;
			}
		}
		log::debug << "Using a " << num_levels << " level heirarchy" << std::endl;
		levels.reserve(num_levels);
		setup_space(num_levels);
		for (auto i : range(num_levels-1)) {
			levels[i].R.associate(&levels[i].P);
			levels[i].P.fine_op = &levels[i].A;
			levels[i].R.set_registry(kreg);
			levels[i].P.set_registry(kreg);

			auto lvl = num_levels - i - 1;
			setup_interp(lvl);
			setup_operator(lvl);
			setup_relax(lvl);
		}
		setup_cg_solve();
	}


	virtual void ncycle(int lvl, grid_func & x, const grid_func & b,
		int n=1)
	{
		auto & A = levels[lvl].A;

		timer_begin("relaxation");
		levels[lvl].presmoother(A, x, b);
		timer_end("relaxation");

		grid_func & residual = levels[lvl].res;
		timer_begin("residual");
		A.residual(x, b, residual);
		timer_end("residual");

		levels[lvl].P.residual = &levels[lvl].res;
		log_residual(lvl, residual);

		auto & coarse_b = levels[lvl+1].b;
		auto & coarse_x = levels[lvl+1].x;
		timer_begin("restrict");
		levels[lvl].R.apply(residual, coarse_b);
		timer_end("restrict");
		coarse_x.set(0.0);

		timer_down();

		int coarse_lvl = levels.size() - 1;
		if (lvl+1 == coarse_lvl) {
			timer_begin("coarse-solve");
			coarse_solver(levels[levels.size()-1].A, coarse_x, coarse_b);
			timer_end("coarse-solve");
		} else {
			for (auto i : range(n)) {
				(void)i;
				ncycle(lvl+1, coarse_x, coarse_b,n);
			}
		}

		timer_up();

		timer_begin("interp-add");
		x += levels[lvl].P * coarse_x;
		timer_end("interp-add");

		timer_begin("relaxation");
		levels[lvl].postsmoother(A, x, b);
		timer_end("relaxation");

		if (log::info.active()) {
			A.residual(x,b,residual);
			log_residual(lvl, residual);
		}
	}


	virtual grid_func solve(const grid_func & b)
	{
		grid_func x = grid_func::zeros_like(b);
		int maxiter = conf->get<int>("solver.max-iter", 10);
		real_t tol = conf->get<real_t>("solver.tol", 1e-8);
		levels[0].A.residual(x,b,levels[0].res);
		real_t res0_l2 = levels[0].res.template lp_norm<2>();
		log::info << "Initial residual l2 norm: " << res0_l2 << std::endl;

		timer_begin("solve");
		for (auto i: range(maxiter)) {
			vcycle(x, b);
			levels[0].A.residual(x,b,levels[0].res);
			real_t res_l2 = levels[0].res.template lp_norm<2>();
			real_t rel_l2 = res_l2 / res0_l2;
			log::status << "Iteration " << i << " relative l2 norm: " << rel_l2 << std::endl;
			if (rel_l2 < tol) break;
		}
		timer_end("solve");

		return x;
	}


	virtual void solve(const grid_func & b, grid_func & x)
	{
		int maxiter = conf->get<int>("solver.max-iter", 10);
		real_t tol = conf->get<real_t>("solver.tol", 1e-8);
		levels[0].A.residual(x,b,levels[0].res);
		real_t res0_l2 = levels[0].res.template lp_norm<2>();
		log::info << "Initial residual l2 norm: " << res0_l2 << std::endl;

		timer_begin("solve");

		for (auto i: range(maxiter)) {
			vcycle(x, b);
			levels[0].A.residual(x,b,levels[0].res);
			real_t res_l2 = levels[0].res.template lp_norm<2>();
			real_t rel_l2 = res_l2 / res0_l2;
			log::status << "Iteration " << i << " relative l2 norm: " << rel_l2 << std::endl;
			if (rel_l2 < tol) break;
		}
		timer_end("solve");
	}

	virtual void vcycle(grid_func & x, const grid_func & b)
	{
		if (levels.size() == 1)
			coarse_solver(levels[0].A, x, b);
		else
			ncycle(0, x, b);
	}

	virtual int compute_num_levels(stencil_op & fop) { return 0; }


	virtual void setup_relax_plane(stencil_op & sop, LevelType & level) {}


	virtual void relax_plane(const stencil_op & so,
	                         grid_func & x,
	                         const grid_func & b,
	                         cycle::Dir cdir,
	                         LevelType & level) {}


	config::reader & get_config() { return *conf; }

protected:
	std::vector<LevelType> levels;
	std::function<void(const discrete_op<grid_func> & A, grid_func &x, const grid_func &b)> coarse_solver;
	std::shared_ptr<config::reader> conf;
	std::shared_ptr<registry> kreg;
	grid_func ABD;
	real_t *bbd;
};

}

#endif
