#ifndef CEDAR_MULTILEVEL_SETTINGS_H
#define CEDAR_MULTILEVEL_SETTINGS_H

#include <iostream>

#include <cedar/types.h>

namespace cedar {

	struct ml_settings
	{
		enum class relax_type {point, line_x, line_y, line_xy,
				plane_xy, plane_xz, plane_yz, plane_xyz};
		enum class cycle_type {v, f};
		static std::map<relax_type, std::string> relax_name;
		friend std::ostream & operator<<(std::ostream & os, const ml_settings & obj);
		void init(config & conf);


		relax_type relaxation;
		cycle_type cycle;
		int nrelax_pre;
		int nrelax_post;
		int num_levels;
		int maxiter;
		real_t tol;
	};

}

#endif
