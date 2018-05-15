#ifndef CEDAR_MULTILEVEL_SETTINGS_H
#define CEDAR_MULTILEVEL_SETTINGS_H

#include <iostream>
#include <vector>

#include <cedar/types.h>
#include <cedar/config.h>

namespace cedar {

	struct redist_settings
	{
		enum class search {astar, coarsen, manual};

		search search_strat;
		std::vector<std::vector<int>> path;
		int min_coarse;
		float machine_latency;
		float machine_bandwidth;
		float machine_fprate;

		friend std::ostream & operator<<(std::ostream & os, const redist_settings & obj);
		void init(config & conf);
	};


	struct ml_settings
	{
		enum class relax_type {point, line_x, line_y, line_xy,
				plane_xy, plane_xz, plane_yz, plane_xyz};
		enum class cycle_type {v, f};
		enum class cg_type {lu, serial, redist};
		static std::map<relax_type, std::string> relax_name;
		friend std::ostream & operator<<(std::ostream & os, const ml_settings & obj);
		void init(config & conf);

		relax_type relaxation;
		cycle_type cycle;
		cg_type coarse_solver;
		std::shared_ptr<config> coarse_config;
		int nrelax_pre;
		int nrelax_post;
		int num_levels;
		int maxiter;
		real_t tol;
		int min_coarse;

		redist_settings rsettings;
	};

}

#endif
