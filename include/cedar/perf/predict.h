#ifndef CEDAR_PERF_PREDICT_H
#define CEDAR_PERF_PREDICT_H

#include <vector>
#include <array>

#include <cedar/types.h>
#include <cedar/perf/perf_factory.h>

namespace cedar {
	template <unsigned short ND>
		std::array<int, ND> choose_redist(config & conf,
		                                  std::array<int, ND> nproc,
		                                  std::array<len_t, ND> nglobal)
	{
		std::array<int, ND> nblocks;
		for (unsigned short i = 0; i < ND; ++i) nblocks[i] = 1;

		auto search_strat = conf.get<std::string>("redist.search.strategy", "astar");

		if (search_strat == "manual") {
			auto path = conf.getnvec<int>("redist.search.path");
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
		} else if (search_strat == "coarsen") {
			for (unsigned short i = 0; i < ND; ++i) {
				nblocks[i] = nproc[i] / 2;
				if (nblocks[i] < 1) nblocks[i] = 1;
			}
		} else if (search_strat == "astar") {
			if (ND == 2) {
				auto model = perf_factory::astar_vcycle(conf, nproc[0], nproc[1], nglobal[0], nglobal[1]);
				nblocks[0] = model->nblocks(0);
				nblocks[1] = model->nblocks(1);
			} else
				log::error << search_strat << " search strategy not implemented for dimension: " << ND << std::endl;
		} else {
			log::error << "Search strategy not implemented: " << search_strat << std::endl;
		}

		return nblocks;
	}
}

#endif
