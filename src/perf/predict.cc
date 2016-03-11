#include <boxmg/perf/perf_factory.h>

#include <boxmg/perf/predict.h>

using namespace boxmg;

std::vector<int> predict_redist(int nprocx, int nprocy,
                                int ngx, int ngy)
{
	std::vector<int> nblocks({1,1});

	auto model = perf_factory::produce_vcycle(nprocx*nprocy, ngx, ngy);

	auto nb = model->get_nchunks();

	for (auto i : range(nb)) {


	}

	return nblocks;
}
