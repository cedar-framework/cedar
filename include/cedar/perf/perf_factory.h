#ifndef CEDAR_PERF_FACTORY_H
#define CEDAR_PERF_FACTORY_H

#include <vector>
#include <memory>

#include <cedar/perf/vcycle_model.h>
#include <cedar/multilevel_settings.h>


namespace cedar {
struct perf_factory
{
	static std::shared_ptr<vcycle_model> produce_vcycle(redist_settings & conf, int npx, int npy, len_t nx, len_t ny);
	static std::shared_ptr<vcycle_model> astar_vcycle(redist_settings & conf, int npx, int npy, len_t nx, len_t ny);
	static std::shared_ptr<vcycle_model> manual_vcycle(redist_settings & conf, int npx, int npy, len_t nx, len_t ny);
	static std::shared_ptr<vcycle_model> random_vcycle(redist_settings & conf, int npx, int npy, len_t nx, len_t ny, std::vector<int> path={});
	static std::shared_ptr<vcycle_model> dfs_vcycle(redist_settings & conf, int npx, int npy, len_t nx, len_t ny, bool terminate=false, int rlevel=0);
	static std::shared_ptr<vcycle_model> produce_vcycle(redist_settings & conf, int npx, int npy, int npz, len_t nx, len_t ny, len_t nz, bool terminate=false);
	static std::array<len_t,2> graph_vcycle(std::ostream & os, int npx, int npy, len_t nx, len_t ny, bool terminate=false, int rlevel = 0);
};
}

#endif
