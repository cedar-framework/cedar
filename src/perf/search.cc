#include <memory>

#include <cedar/perf/redist_generator.h>
#include <cedar/perf/greedy_iterator.h>
#include <cedar/perf/perf_factory.h>
#include <cedar/perf/const_model.h>
#include <cedar/perf/search.h>

using namespace cedar;

bool perf_problem::goal_test(perf_state & state)
{
	return (state.model->grid(0).nproc() == 1);
}


perf_state perf_problem::result(perf_state & state, std::array<int,2> action)
{
	perf_state ret;

	auto model = state.model;
	model->nblocks(0) = action[0];
	model->nblocks(1) = action[1];
	auto & topoc = model->grid(0);

	auto cg_model = perf_factory::produce_vcycle(settings, action[0], action[1],
	                                             topoc.nglobal(0), topoc.nglobal(1));

	model->set_cgperf(cg_model);
	ret.model = cg_model;

	return ret;
}


float perf_problem::step_cost(perf_state & state, std::array<int,2> /*act*/)
{
	auto model = state.model;

	return model->tcgsolve();
}


std::vector<std::array<int,2>> perf_problem::actions(perf_state & state)
{
	std::vector<std::array<int,2>> ret;

	auto model = state.model;
	auto & topoc = model->grid(0);
	auto npx = topoc.nproc(0);
	auto npy = topoc.nproc(1);

	int min_coarse = settings.min_coarse;

	auto redist_subsets = redist_generator<greedy_iterator>({npx, npy}, {topoc.nglobal(0), topoc.nglobal(1)}, min_coarse);
	for (auto nblocks : redist_subsets) {
		ret.push_back(nblocks);
	}

	return ret;
}


namespace cedar {
	bool operator<(const perf_state & s1, const perf_state & s2)
	{
		// TODO: handle this recursively
		// if (s1.model == nullptr) return true;
		// else if (s2.model == nullptr) return false;

		len_t ngc[2][2];
		len_t np[2][2];

		for (auto i = 0; i < 2; i++) {
			ngc[0][i] = s1.model->grid(0).nglobal(i);
			ngc[1][i] = s2.model->grid(0).nglobal(i);
			np[0][i] = s1.model->grid(0).nproc(i);
			np[1][i] = s2.model->grid(0).nproc(i);
		}

		if (ngc[0][0] < ngc[1][0]) return true;
		else if (ngc[0][0] == ngc[1][0] and ngc[0][1] < ngc[1][1]) return true;
		else if (ngc[0][0] == ngc[1][0] and ngc[0][1] == ngc[1][1] and np[0][0] < np[1][0]) return true;
		else if (ngc[0][0] == ngc[1][0] and ngc[0][1] == ngc[1][1] and np[0][0] == np[1][0] and np[0][1] < np[1][1]) return true;
		else return false;
	}
}
