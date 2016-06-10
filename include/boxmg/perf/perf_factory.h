#ifndef BOXMG_PERF_FACTORY_H
#define BOXMG_PERF_FACTORY_H

#include <vector>
#include <memory>

#include <boxmg/perf/vcycle_model.h>


namespace boxmg {
struct perf_factory
{
	static std::shared_ptr<vcycle_model> produce_vcycle(config::reader & conf, int npx, int npy, len_t nx, len_t ny);
	static std::shared_ptr<vcycle_model> astar_vcycle(config::reader & conf, int npx, int npy, len_t nx, len_t ny);
	static std::shared_ptr<vcycle_model> manual_vcycle(config::reader & conf, int npx, int npy, len_t nx, len_t ny);
	static std::shared_ptr<vcycle_model> random_vcycle(config::reader & conf, int npx, int npy, len_t nx, len_t ny);
	static std::shared_ptr<vcycle_model> dfs_vcycle(config::reader & conf, int npx, int npy, len_t nx, len_t ny, bool terminate=false, int rlevel=0);
	static std::shared_ptr<vcycle_model> produce_vcycle(config::reader & conf, int npx, int npy, int npz, len_t nx, len_t ny, len_t nz, bool terminate=false);
	static std::array<len_t,2> graph_vcycle(std::ostream & os, int npx, int npy, len_t nx, len_t ny, bool terminate=false, int rlevel = 0);
};
}

#endif
