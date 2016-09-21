#ifndef BOXMG_PERF_PREDICT_H
#define BOXMG_PERF_PREDICT_H

#include <vector>

#include <boxmg/types.h>

namespace boxmg {
	std::vector<int> predict_redist(config::reader & conf,
	                                int nprocx, int nprocy,
	                                len_t ngx, len_t ngy);
}

#endif
