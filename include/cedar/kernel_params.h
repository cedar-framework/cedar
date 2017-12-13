#ifndef CEDAR_KERNEL_PARAMS_H
#define CEDAR_KERNEL_PARAMS_H

#include <array>
#include <memory>

#include <cedar/config/reader.h>

namespace cedar
{
	struct kernel_params
	{
		std::array<bool,3> periodic;
		bool relax_symmetric;
		bool definite;
		bool ml_relax;

		int per_mask() const {
			int mask = 0;
			for (int i = 0; i < 3; ++i) {
				if (periodic[i])
					mask |= (1<<i);
			}
			return mask;
		}
	};


	std::shared_ptr<kernel_params> build_kernel_params(config::reader & conf);
}

#endif
