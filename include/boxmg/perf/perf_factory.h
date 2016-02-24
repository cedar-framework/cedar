#ifndef BOXMG_PERF_FACTORY_H
#define BOXMG_PERF_FACTORY_H

#include <vector>
#include <memory>

#include <boxmg/perf/vcycle_model.h>


namespace boxmg {
struct perf_factory
{
	static std::shared_ptr<vcycle_model> produce_vcycle(int np, len_t nx, len_t ny, bool terminate=false);
	static std::shared_ptr<vcycle_model> produce_vcycle(int np, len_t nx, len_t ny, len_t nz, bool terminate=false);
};
}

#endif
