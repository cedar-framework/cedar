#ifndef CEDAR_MULTILEVEL_H
#define CEDAR_MULTILEVEL_H

#include <cedar/types.h>
#include <cedar/kernel_registry.h>
#include <cedar/util/timer.h>
#include <cedar/cycle/types.h>
#include <cedar/cycle/vcycle.h>
#include <cedar/cycle/fcycle.h>

namespace cedar {


	/**
	   Base class for multilevel solvers.

	   Base class for every mutilevel solve class.  This class
	   performs operations generically by taking a structured listing
	   of types and kernel implementations (registry) as input.

	   @tparam level_container Generic container for data that is stored on every level.
	   @tparam fsten Stencil used for the fine-grid operator.
	   @tparam child Class that inherits this.
	*/
template <exec_mode emode, class level_container, class fsten, class child>
class multilevel
{
public:
	template<class sten>
		using level_t = typename level_container::template level_t<sten>;
	template<class sten>
		using stencil_op = typename level_t<fsten>::template stencil_op<sten>;
	using grid_func = typename level_t<fsten>::grid_func;
	// extract kernel names for convenience
	using stypes = typename level_t<fsten>::stypes;
	template<relax_dir rdir>
		using line_relax = kernels::line_relax<stypes, rdir>;
	template<relax_dir rdir>
		using plane_relax = kernels::plane_relax<stypes, rdir>;
	using point_relax = kernels::point_relax<stypes>;
	using coarsen_op = kernels::coarsen_op<stypes>;
	using interp_add = kernels::interp_add<stypes>;
	using matvec = kernels::matvec<stypes>;
	using residual = kernels::residual<stypes>;
	using setup_prolong = kernels::setup_interp<stypes>;
	using solve_cg = kernels::solve_cg<stypes>;
	using kern_manager = kernel_manager<klist<stypes, emode>>;

	using conf_ptr = std::shared_ptr<config::reader>;
multilevel(stencil_op<fsten> & fop) : levels(fop), conf(std::make_shared<config::reader>("config.json")) {
		setup_cycler();
	}


multilevel(stencil_op<fsten> & fop, conf_ptr cfg): levels(fop), conf(cfg) {
		setup_cycler();
	}


	virtual ~multilevel() {}


	void setup_cycler()
	{
		auto cycle_type = conf->get<std::string>("solver.cycle.type", "v");
		if (cycle_type == "v") {
			this->cycle = std::make_unique<cycle::vcycle<emode, level_container, fsten>>(
				levels, coarse_solver);
		} else if (cycle_type == "f") {
			this->cycle = std::make_unique<cycle::fcycle<emode, level_container, fsten>>(
				levels, coarse_solver);
		} else {
			log::error << "Cycle type: " << cycle_type << " not recognized!" << std::endl;
		}
	}


	void vcycle(grid_func & x, const grid_func & b)
	{
		auto cycle_type = conf->get<std::string>("solver.cycle.type", "v");
		assert(cycle_type == "v");
		this->cycle->run(x, b);
	}


	std::shared_ptr<kern_manager> get_kernels()
	{
		return kman;
	}

	std::size_t nlevels() { return levels.size(); }


	virtual void setup_cg_solve()
	{
		auto & cop = levels.get(levels.size() - 1).A;
		kman->template setup<solve_cg>(cop, ABD);
		auto kernels = get_kernels();
		coarse_solver = [&, kernels](grid_func &x, const grid_func & b) {
			kernels->template run<solve_cg>(x, b, ABD, bbd);
		};
	}


	void setup_interp(std::size_t lvl)
	{
		auto & P = levels.get(lvl+1).P;
		auto & cop = levels.get(lvl+1).A;
		if (lvl == 0) {
			auto & fop = levels.template get<fsten>(lvl).A;
			kman->template run<setup_prolong>(fop, cop, P);
		} else {
			auto & fop = levels.get(lvl).A;
			kman->template run<setup_prolong>(fop, cop, P);
		}
	}


	void setup_operator(std::size_t lvl)
	{
		auto & P = levels.get(lvl+1).P;
		auto & cop = levels.get(lvl+1).A;

		if (lvl == 0) {
			auto & fop = levels.template get<fsten>(lvl).A;
			kman->template run<coarsen_op>(P, fop, cop);
		} else {
			auto & fop = levels.get(lvl).A;
			kman->template run<coarsen_op>(P, fop, cop);
		}
	}

	template<class sten>
		void setup_relax_helper(level_t<sten> & level, std::size_t lvl)
	{
		auto & sop = level.A;

		std::string relax_type = conf->get<std::string>("solver.relaxation", "point");
		int nrelax_pre = conf->get<int>("solver.cycle.nrelax-pre", 2);
		int nrelax_post = conf->get<int>("solver.cycle.nrelax-post", 1);

		if (relax_type == "point")
			kman->template setup<point_relax>(sop, level.SOR[0]);
		else if (relax_type == "line-x")
			kman->template setup<line_relax<relax_dir::x>>(sop, level.SOR[0]);
		else if (relax_type == "line-y")
			kman->template setup<line_relax<relax_dir::y>>(sop, level.SOR[0]);
		else if (relax_type == "line-xy") {
			kman->template setup<line_relax<relax_dir::x>>(sop, level.SOR[0]);
			kman->template setup<line_relax<relax_dir::y>>(sop, level.SOR[1]);
		}
		else if (relax_type == "plane") {
			kman->template setup<plane_relax<relax_dir::xy>>(sop);
			kman->template setup<plane_relax<relax_dir::xz>>(sop);
			kman->template setup<plane_relax<relax_dir::yz>>(sop);
		}
		else if (relax_type == "plane-xy")
			kman->template setup<plane_relax<relax_dir::xy>>(sop);
		else if (relax_type == "plane-xz")
			kman->template setup<plane_relax<relax_dir::xz>>(sop);
		else if (relax_type == "plane-yz")
			kman->template setup<plane_relax<relax_dir::yz>>(sop);
		else
			log::error << "Invalid relaxation: " << relax_type << std::endl;

		auto kernels = get_kernels();

		level.presmoother = [&,lvl,nrelax_pre,kernels,relax_type](const stencil_op<sten> &A,
		                                                          grid_func &x, const grid_func&b) {
			for (auto i : range(nrelax_pre)) {
				(void) i;
				if (relax_type == "point")
					kernels->template run<point_relax>(A, x, b, level.SOR[0], cycle::Dir::DOWN);
				else if (relax_type == "line-x")
					kernels->template run<line_relax<relax_dir::x>>(A, x, b, level.SOR[0], level.res, cycle::Dir::DOWN);
				else if (relax_type == "line-y")
					kernels->template run<line_relax<relax_dir::y>>(A, x, b, level.SOR[0], level.res, cycle::Dir::DOWN);
				else if (relax_type == "line-xy") {
					kernels->template run<line_relax<relax_dir::x>>(A, x, b, level.SOR[0], level.res, cycle::Dir::DOWN);
					kernels->template run<line_relax<relax_dir::y>>(A, x, b, level.SOR[1], level.res, cycle::Dir::DOWN);
				}
				else if (relax_type == "plane") {
					kernels->template run<plane_relax<relax_dir::xy>>(A, x, b, cycle::Dir::DOWN);
					kernels->template run<plane_relax<relax_dir::yz>>(A, x, b, cycle::Dir::DOWN);
					kernels->template run<plane_relax<relax_dir::xz>>(A, x, b, cycle::Dir::DOWN);
				}
				else if (relax_type == "plane-xy")
					kernels->template run<plane_relax<relax_dir::xy>>(A, x, b, cycle::Dir::DOWN);
				else if (relax_type == "plane-xz")
					kernels->template run<plane_relax<relax_dir::xz>>(A, x, b, cycle::Dir::DOWN);
				else if (relax_type == "plane-yz")
					kernels->template run<plane_relax<relax_dir::yz>>(A, x, b, cycle::Dir::DOWN);
				else
					log::error << "Invalid relaxation: " << relax_type << std::endl;
			}
		};
		level.postsmoother = [&,lvl,nrelax_post,kernels,relax_type](const stencil_op<sten> &A,
		                                                            grid_func &x, const grid_func&b) {
			for (auto i: range(nrelax_post)) {
				(void) i;
				if (relax_type == "point")
					kernels->template run<point_relax>(A, x, b, level.SOR[0], cycle::Dir::UP);
				else if (relax_type == "line-x")
					kernels->template run<line_relax<relax_dir::x>>(A, x, b, level.SOR[0], level.res, cycle::Dir::UP);
				else if (relax_type == "line-y")
					kernels->template run<line_relax<relax_dir::y>>(A, x, b, level.SOR[0], level.res, cycle::Dir::UP);
				else if (relax_type == "line-xy") {
					kernels->template run<line_relax<relax_dir::y>>(A, x, b, level.SOR[1], level.res, cycle::Dir::UP);
					kernels->template run<line_relax<relax_dir::x>>(A, x, b, level.SOR[0], level.res, cycle::Dir::UP);
				}
				else if (relax_type == "plane") {
					kernels->template run<plane_relax<relax_dir::xz>>(A, x, b, cycle::Dir::UP);
					kernels->template run<plane_relax<relax_dir::yz>>(A, x, b, cycle::Dir::UP);
					kernels->template run<plane_relax<relax_dir::xy>>(A, x, b, cycle::Dir::UP);
				}
				else if (relax_type == "plane-xy")
					kernels->template run<plane_relax<relax_dir::xy>>(A, x, b, cycle::Dir::UP);
				else if (relax_type == "plane-xz")
					kernels->template run<plane_relax<relax_dir::xz>>(A, x, b, cycle::Dir::UP);
				else if (relax_type == "plane-yz")
					kernels->template run<plane_relax<relax_dir::yz>>(A, x, b, cycle::Dir::UP);
				else
					log::error << "Invalid relaxation: " << relax_type << std::endl;
			}
		};
	}

	void setup_relax(std::size_t lvl)
	{
		if (lvl == 0) {
			auto & level = levels.template get<fsten>(lvl);
			setup_relax_helper(level, lvl);
		} else {
			auto & level = levels.get(lvl);
			setup_relax_helper(level, lvl);
		}
	}


	void setup_space(std::size_t nlevels)
	{
		static_cast<child*>(this)->setup_space(nlevels);
	}


	void setup(stencil_op<fsten> & fop)
	{
		this->cycle->set_kernels(this->kman);
		auto num_levels = compute_num_levels(fop);
		auto nlevels_conf = conf->get<int>("solver.num-levels", -1);
		if (nlevels_conf > 0) {
			if (static_cast<std::size_t>(nlevels_conf) > num_levels) {
				log::error << "too many levels specified" << std::endl;
			} else {
				num_levels = nlevels_conf;
			}
		}
		log::debug << "Using a " << num_levels << " level heirarchy" << std::endl;
		setup_space(num_levels);
		timer_begin("setup");
		for (std::size_t i = 0; i < num_levels - 1; ++i) {
			setup_interp(i);
			setup_operator(i);
			setup_relax(i);
		}
		setup_cg_solve();
		timer_end("setup");
	}


	virtual grid_func solve(const grid_func & b)
	{
		grid_func x = grid_func::zeros_like(b);

		solve(b, x);

		return x;
	}


	virtual void solve(const grid_func & b, grid_func & x)
	{
		auto & level = levels.template get<fsten>(0);
		int maxiter = conf->get<int>("solver.max-iter", 10);
		real_t tol = conf->get<real_t>("solver.tol", 1e-8);
		kman->template run<residual>(level.A,x,b,level.res);
		real_t res0_l2 = level.res.template lp_norm<2>();
		log::info << "Initial residual l2 norm: " << res0_l2 << std::endl;

		timer_begin("solve");

		for (auto i: range(maxiter)) {
			cycle->run(x, b);
			kman->template run<residual>(level.A,x,b,level.res);
			real_t res_l2 = level.res.template lp_norm<2>();
			real_t rel_l2 = res_l2 / res0_l2;
			log::status << "Iteration " << i << " relative l2 norm: " << rel_l2 << std::endl;
			if (rel_l2 < tol) break;
		}
		timer_end("solve");
	}


	std::size_t compute_num_levels(stencil_op<fsten> & fop)
	{
		return static_cast<child*>(this)->compute_num_levels(fop);
	}


	config::reader & get_config() { return *conf; }
	level_container levels;

protected:
	std::unique_ptr<cycle::cycle<emode, level_container, fsten>> cycle;
	std::function<void(grid_func &x, const grid_func &b)> coarse_solver;
	std::shared_ptr<config::reader> conf;
	std::shared_ptr<kern_manager> kman;
	grid_func ABD;
	real_t *bbd;
};

}

#endif
