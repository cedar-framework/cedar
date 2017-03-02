#ifndef BOXMG_PERF_PREDICT_H
#define BOXMG_PERF_PREDICT_H

#include <vector>
#include <array>

#include <boxmg/types.h>

namespace boxmg {
	std::vector<int> predict_redist(config::reader & conf,
	                                int nprocx, int nprocy,
	                                len_t ngx, len_t ngy);

	template <unsigned short ND>
		std::array<int, ND> choose_redist(config::reader & conf,
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
		} else {
			log::error << "Search strategy not implemented: " << search_strat << std::endl;
		}

		return nblocks;
	}
}

#endif
