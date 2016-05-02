#include <boxmg/perf/perf_factory.h>

#include <boxmg/perf/predict.h>

namespace boxmg {

std::vector<int> predict_redist(int nprocx, int nprocy,
                                len_t ngx, len_t ngy, int nlevels)
{
	std::vector<int> nblocks({1,1});

	auto model = perf_factory::produce_vcycle(nprocx, nprocy, ngx, ngy, false, nlevels);

	nblocks[0] = model->nblocks(0);
	nblocks[1] = model->nblocks(1);

	log::status << '\n' << *model << std::endl;

	return nblocks;
}


}
