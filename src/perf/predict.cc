#include <boxmg/perf/perf_factory.h>

#include <boxmg/perf/predict.h>

using namespace boxmg;

std::vector<int> predict_redist(int nprocx, int nprocy,
                                int ngx, int ngy)
{
	std::vector<int> nblocks({1,1});

	auto model = perf_factory::produce_vcycle(nprocx*nprocy, ngx, ngy);

	auto nb = model->get_nchunks();

	// assuming number of blocks is a power of 2
	int nref = std::log2(nb);

	// basing on nproc for now (instead of nglobal)
	// greedy assignment of processors to blocks
	for (auto i : range(nref)) {
		if ((nprocx / nblocks[0]) > (nprocy / nblocks[1]))
			nblocks[0] *= 2;
		else
			nblocks[1] *= 2;
	}

	return nblocks;
}
