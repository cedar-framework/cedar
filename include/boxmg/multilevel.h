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
multilevel() : conf("config.json") {};
	virtual ~multilevel() {delete[] bbd;}

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
		kreg->setup_interp(lvl, lvl-1, levels.size(), fop, cop, P);
	}


	virtual void setup_operator(int lvl)
	{
		auto & P = level(lvl).P;
		auto & fop = level(lvl).A;
		auto & cop = level(lvl-1).A;
		kreg->galerkin_prod(lvl, lvl-1, levels.size(), P, fop, cop);
	}


	virtual void setup_relax(int lvl)
	{
		auto & sop = level(lvl).A;

		std::string relax_type = conf.get<std::string>("solver.relaxation", "point");
		int nrelax_pre = conf.get<int>("solver.cycle.nrelax-pre", 2);
		int nrelax_post = conf.get<int>("solver.cycle.nrelax-post", 1);
		auto kernels = kernel_registry();

		if (relax_type == "point")
			kernels->setup_relax(sop, level(lvl).SOR[0]);
		else if (relax_type == "line-x")
			kernels->setup_relax_x(sop, level(lvl).SOR[0]);
		else if (relax_type == "line-y")
			kernels->setup_relax_y(sop, level(lvl).SOR[0]);
		else { // line-xy
			kernels->setup_relax_x(sop, level(lvl).SOR[0]);
			kernels->setup_relax_y(sop, level(lvl).SOR[1]);
		}

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
				else {
					kernels->relax_lines_x(av, x, b, level(lvl).SOR[0], level(lvl).res, cycle::Dir::DOWN);
					kernels->relax_lines_y(av, x, b, level(lvl).SOR[1], level(lvl).res, cycle::Dir::DOWN);
				}
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
				else {
					kernels->relax_lines_y(av, x, b, level(lvl).SOR[1], level(lvl).res, cycle::Dir::UP);
					kernels->relax_lines_x(av, x, b, level(lvl).SOR[0], level(lvl).res, cycle::Dir::UP);
				}
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

		levels[lvl].presmoother(A, x, b);

		grid_func & residual = levels[lvl].res;
		A.residual(x, b, residual);

		levels[lvl].P.residual = &levels[lvl].res;
		log_residual(lvl, residual);

		auto & coarse_b = levels[lvl+1].b;
		auto & coarse_x = levels[lvl+1].x;
		levels[lvl].R.apply(residual, coarse_b);
		coarse_x.set(0.0);

		if (lvl == levels.size() - 2) {
			coarse_solver(levels[levels.size()-1].A, coarse_x, coarse_b);
		} else {
			for (auto i : range(n)) {
				(void)i;
				ncycle(lvl+1, coarse_x, coarse_b,n);
			}
		}

		x += levels[lvl].P * coarse_x;

		levels[lvl].postsmoother(A, x, b);

		if (log::info.active()) {
			A.residual(x,b,residual);
			log_residual(lvl, residual);
		}
	}


	virtual grid_func solve(const grid_func & b)
	{
		grid_func x = grid_func::zeros_like(b);
		int maxiter = config::get<int>("solver.max-iter", 10);
		real_t tol = config::get<real_t>("solver.tol", 1e-8);
		levels[0].A.residual(x,b,levels[0].res);
		real_t res0_l2 = levels[0].res.template lp_norm<2>();
		log::info << "Initial residual l2 norm: " << res0_l2 << std::endl;

		timer solve_timer("Solve");
		solve_timer.begin();

		for (auto i: range(maxiter)) {
			ncycle(0, x, b);
			levels[0].A.residual(x,b,levels[0].res);
			real_t res_l2 = levels[0].res.template lp_norm<2>();
			real_t rel_l2 = res_l2 / res0_l2;
			log::status << "Iteration " << i << " relative l2 norm: " << rel_l2 << std::endl;
			if (rel_l2 < tol) break;
		}

		solve_timer.end();

		return x;
	}


	virtual void solve(const grid_func & b, grid_func & x)
	{
		int maxiter = config::get<int>("solver.max-iter", 10);
		real_t tol = config::get<real_t>("solver.tol", 1e-8);
		levels[0].A.residual(x,b,levels[0].res);
		real_t res0_l2 = levels[0].res.template lp_norm<2>();
		log::info << "Initial residual l2 norm: " << res0_l2 << std::endl;

		timer solve_timer("Solve");
		solve_timer.begin();

		for (auto i: range(maxiter)) {
			ncycle(0, x, b);
			levels[0].A.residual(x,b,levels[0].res);
			real_t res_l2 = levels[0].res.template lp_norm<2>();
			real_t rel_l2 = res_l2 / res0_l2;
			log::status << "Iteration " << i << " relative l2 norm: " << rel_l2 << std::endl;
			if (rel_l2 < tol) break;
		}

		solve_timer.end();
	}


	virtual int compute_num_levels(stencil_op & fop) { return 0; }

protected:
	std::vector<LevelType> levels;
	std::function<void(const discrete_op<grid_func> & A, grid_func &x, const grid_func &b)> coarse_solver;
	config::Reader conf;
	std::shared_ptr<registry> kreg;
	grid_func ABD;
	real_t *bbd;
};

}

#endif
