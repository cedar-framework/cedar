#ifndef CEDAR_KERNEL_PARAMS_H
#define CEDAR_KERNEL_PARAMS_H

#include <array>
#include <memory>

#include <cedar/config/reader.h>

namespace cedar
{
	struct ml_relax_params
	{
		ml_relax_params() : enabled(false), shm(false), min_gsz(3) {}
		void init(config::reader & conf);

		bool enabled; /** Whether multilevel line relaxation is enabled */
		bool shm;     /** Whether shared memory optimization is enabled */
		int min_gsz;  /** Coarsening factor for multilevel line relaxation */
	};


	struct kernel_params
	{
		std::array<bool,3> periodic;
		bool relax_symmetric;
		bool definite;
		ml_relax_params ml_relax;

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
