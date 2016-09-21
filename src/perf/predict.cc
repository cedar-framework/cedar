#include <boxmg/perf/perf_factory.h>

#include <boxmg/perf/predict.h>

namespace boxmg {

std::vector<int> predict_redist(config::reader & conf, int nprocx, int nprocy,
                                len_t ngx, len_t ngy)
{
	std::vector<int> nblocks({1,1});

	auto search_strat = conf.get<std::string>("redist.search.strategy", "astar");

	if (search_strat == "astar") {
		auto model = perf_factory::astar_vcycle(conf, nprocx, nprocy, ngx, ngy);
		nblocks[0] = model->nblocks(0);
		nblocks[1] = model->nblocks(1);
		log::status << '\n' << *model << std::endl;
	} else if (search_strat == "manual") {
		auto path = conf.getnvec<int>("redist.search.path");
		bool found = false;
		for (auto & proc : path) {
			if (found) {
				nblocks[0] = proc[0];
				nblocks[1] = proc[1];
				break;
			} else if (proc[0] == nprocx and proc[1] == nprocy) {
				found = true;
			}
		}
	} else {
		auto model = perf_factory::dfs_vcycle(conf, nprocx, nprocy, ngx, ngy, false);
		nblocks[0] = model->nblocks(0);
		nblocks[1] = model->nblocks(1);
		log::status << '\n' << *model << std::endl;
	}

	return nblocks;
}


}
