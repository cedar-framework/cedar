#ifndef CEDAR_PERF_PREDICT_H
#define CEDAR_PERF_PREDICT_H

#include <vector>
#include <array>

#include <cedar/types.h>
#include <cedar/perf/perf_factory.h>

namespace cedar {
	template <unsigned short ND>
		std::array<int, ND> choose_redist(redist_settings & settings,
		                                  std::array<int, ND> nproc,
		                                  std::array<len_t, ND> nglobal)
	{
		std::array<int, ND> nblocks;
		for (unsigned short i = 0; i < ND; ++i) nblocks[i] = 1;

		if (settings.search_strat == redist_settings::search::manual) {
			auto & path = settings.path;
			bool found = false;
			for (auto & proc : path) {
				if (found) {
					for (unsigned short i = 0; i < ND; ++i)
						nblocks[i] = proc[i];
					break;
				} else {
					bool cond = true;
					for (unsigned short i = 0; i < ND; ++i)
						cond = cond and (nproc[i] == proc[i]);
					if (cond)
						found = true;
				}
			}
		} else if (settings.search_strat == redist_settings::search::coarsen) {
			float curr_min = std::numeric_limits<float>::max();
			int curr_ind = 0;
			for (unsigned short i = 0; i < ND; ++i) {
				nblocks[i] = nproc[i];
				float val = nglobal[i] / nproc[i];
				if (nproc[i] > 1 and val < curr_min) {
					curr_min = val;
					curr_ind = i;
				}
			}
			nblocks[curr_ind] = nblocks[curr_ind] / 2;
			if (nblocks[curr_ind] < 1) nblocks[curr_ind] = 1;
		} else if (settings.search_strat == redist_settings::search::astar) {
			if (ND == 2) {
				auto model = perf_factory::astar_vcycle(settings, nproc[0], nproc[1], nglobal[0], nglobal[1]);
				nblocks[0] = model->nblocks(0);
				nblocks[1] = model->nblocks(1);
			} else
				log::error << "astar search strategy not implemented for dimension: " << ND << std::endl;
		} else {
			log::error << "Search strategy not implemented" << std::endl;
		}

		return nblocks;
	}
}

#endif
