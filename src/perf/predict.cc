#include <boxmg/perf/perf_factory.h>

#include <boxmg/perf/predict.h>

namespace boxmg {

std::vector<int> predict_redist(int nprocx, int nprocy,
                                int ngx, int ngy)
{
	std::vector<int> nblocks({1,1});

	auto model = perf_factory::produce_vcycle(nprocx, nprocy, ngx, ngy);

	nblocks[0] = model->nblocks(0);
	nblocks[1] = model->nblocks(1);

	return nblocks;
}


}
