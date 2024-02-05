#ifndef CEDAR_KERNEL_PARAMS_H
#define CEDAR_KERNEL_PARAMS_H

#include <array>
#include <memory>

#include <cedar/config.h>

namespace cedar
{
	struct ml_relax_params
	{
		ml_relax_params() : enabled(false), min_gsz(3) {}
		void init(config & conf);

		bool enabled; /** Whether multilevel line relaxation is enabled */
		int min_gsz;  /** Coarsening factor for multilevel line relaxation */
		bool factorize; /** Factorize local blocks (compared to elimintation */
	};


	struct kernel_params
	{
		enum class halo_lib {tausch, msg};
		friend std::ostream & operator<<(std::ostream & os, const kernel_params & obj);

		std::array<bool,3> periodic;
		bool relax_symmetric;
		bool definite;
		halo_lib halo;
		ml_relax_params ml_relax;
		std::shared_ptr<config> plane_config;
		bool plane_agg;
            bool use_gpu;
            bool use_gpu_cholesky;

		int per_mask() const {
			int mask = 0;
			for (int i = 0; i < 3; ++i) {
				if (periodic[i])
					mask |= (1<<i);
			}
			return mask;
		}
	};


	std::shared_ptr<kernel_params> build_kernel_params(config & conf);
}

#endif
