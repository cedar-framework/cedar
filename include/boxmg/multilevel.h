#ifndef BOXMG_MULTILEVEL_H
#define BOXMG_MULTILEVEL_H

#include "boxmg/types.h"
#include "boxmg/util/timer.h"
#include "boxmg/discrete_op.h"

namespace boxmg {

template <class LevelType,class grid_func,class registry>
class multilevel
{

public:
multilevel() : conf("config.json") {};
	virtual ~multilevel() {};

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

protected:
	std::vector<LevelType> levels;
	std::function<void(const discrete_op<grid_func> & A, grid_func &x, const grid_func &b)> coarse_solver;
	config::Reader conf;
	std::shared_ptr<registry> kreg;
};

}

#endif
