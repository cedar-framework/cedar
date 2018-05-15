#include <cedar/multilevel_settings.h>

namespace cedar {

	std::map<ml_settings::relax_type,
	         std::string> ml_settings::relax_name = {{ml_settings::relax_type::point, "point"},
	                                                 {ml_settings::relax_type::line_x, "line-x"},
	                                                 {ml_settings::relax_type::line_y, "line-y"},
	                                                 {ml_settings::relax_type::line_xy, "line-xy"},
	                                                 {ml_settings::relax_type::plane_xy, "plane-xy"},
	                                                 {ml_settings::relax_type::plane_xz, "plane-xz"},
	                                                 {ml_settings::relax_type::plane_yz, "plane-yz"},
	                                                 {ml_settings::relax_type::plane_xyz, "plane-xyz"}};

	void ml_settings::init(config & conf)
	{
		auto relax_type = conf.get<std::string>("solver.relaxation", "point");
		{
			bool found = false;
			for (auto & rmap : relax_name) {
				if (relax_type == rmap.second) {
					found = true;
					this->relaxation = rmap.first;
				}
			}
			if (not found)
				log::error << "invalid relaxation type: " << relax_type << std::endl;
		}

		auto cycle_name = conf.get<std::string>("solver.cycle.type", "v");
		if (cycle_name == "v")
			this->cycle = cycle_type::v;
		else if (cycle_name == "f")
			this->cycle = cycle_type::f;
		else
			log::error << "invalid cycle type: " << cycle_name << std::endl;

		this->nrelax_pre = conf.get<int>("solver.cycle.nrelax-pre", 2);
		this->nrelax_post = conf.get<int>("solver.cycle.nrelax-post", 1);
		this->num_levels = conf.get<int>("solver.num-levels", -1);
		this->maxiter = conf.get<int>("solver.max-iter", 10);
		this->tol = conf.get<real_t>("solver.tol", 1e-8);
		this->min_coarse = conf.get<int>("solver.min_coarse", 3);

		auto cg_solver_str = conf.get<std::string>("solver.cg-solver", "LU");
		if (cg_solver_str == "LU")
			this->coarse_solver = cg_type::lu;
		else if (cg_solver_str == "cedar")
			this->coarse_solver = cg_type::serial;
		else if (cg_solver_str == "redist")
			this->coarse_solver = cg_type::redist;
		else
			log::error << "invalid value for solver.cg-solver" << std::endl;

		this->coarse_config = conf.getconf("cg-config");
		if (!this->coarse_config)
			this->coarse_config = conf.getconf("");

		if (this->coarse_solver == cg_type::redist)
			this->rsettings.init(conf);
	}


	std::ostream & operator<<(std::ostream & os, const ml_settings & obj)
	{
		os << '\n';
		os << "-------------------\n";
		os << "Multilevel Settings\n";
		os << "-------------------\n";

		os << "coarse solver:   ";
		if (obj.coarse_solver == ml_settings::cg_type::lu)
			os << "cholesky\n";
		else if (obj.coarse_solver == ml_settings::cg_type::redist)
			os << "redist\n";
		else
			os << "serial cedar\n";

		os << "relaxation:      " << ml_settings::relax_name[obj.relaxation] << '\n';
		os << "cycle:           ";
		if (obj.cycle == ml_settings::cycle_type::v)
			os << "V";
		else
			os << "F";
		os << '\n';

		os << "min coarse:      " << obj.min_coarse << '\n';
		if (obj.num_levels > 0)
			os << "num levels:     " << obj.num_levels << '\n';
		os << "nrelax pre:      " << obj.nrelax_pre << '\n';
		os << "nrelax post:     " << obj.nrelax_post << '\n';
		os << "maxiter:         " << obj.maxiter << '\n';
		os << "tol:             " << obj.tol;

		if (obj.coarse_solver == ml_settings::cg_type::redist)
			os << obj.rsettings;

		return os;
	}


	void redist_settings::init(config & conf)
	{
		auto search_strat_str = conf.get<std::string>("redist.search.strategy", "coarsen");
		if (search_strat_str == "manual")
			this->search_strat = search::manual;
		else if (search_strat_str == "coarsen")
			this->search_strat = search::coarsen;
		else if (search_strat_str == "astar")
			this->search_strat = search::astar;
		else
			log::error << "Search strategy not recognized: " << search_strat_str << std::endl;


		if (search_strat == search::manual) {
			this->path = conf.getnvec<int>("redist.search.path");
		}

		if (search_strat == search::astar) {
			this->min_coarse = conf.get<int>("solver.min-coarse");
			this->machine_bandwidth = conf.get<float>("machine.bandwidth");
			this->machine_latency = conf.get<float>("machine.latency");
			this->machine_fprate = conf.get<float>("machine.fp_perf");
		}
	}


	std::ostream & operator<<(std::ostream & os, const redist_settings & obj)
	{
		os << '\n';
		os << "redist strategy: ";
		if (obj.search_strat == redist_settings::search::manual)
			os << "manual\n";
		else if (obj.search_strat == redist_settings::search::coarsen)
			os << "coarsen\n";
		else if (obj.search_strat == redist_settings::search::astar)
			os << "optimized\n";

		return os;
	}
}
