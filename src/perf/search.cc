#include <memory>

#include <boxmg/perf/perf_factory.h>
#include <boxmg/perf/const_model.h>
#include <boxmg/perf/search.h>

using namespace boxmg;

static bool keep_refining(int npx, int npy, int nbx, int nby, len_t nlx, len_t nly,
                          int min_coarse)
{
	bool ret = ((npx / nbx) > 0 and (npy / nby) > 0);
	ret = ret and ((npx / nbx) * (npy / nby) > 1);
	ret = ret and (nlx > 2*min_coarse);
	ret = ret and (nly > 2*min_coarse);

	return ret;
}


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

	auto cg_model = perf_factory::produce_vcycle(action[0], action[1],
	                                             topoc.nglobal(0), topoc.nglobal(1));

	// set coarse solve to 0
	model->set_cgperf(cg_model);
	cg_model->set_cgperf(std::make_shared<const_model>(0));
	ret.model = cg_model;

	return ret;
}


float perf_problem::step_cost(perf_state & state, std::array<int,2> act)
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
	len_t nlx = topoc.nglobal(0);
	len_t nly = topoc.nglobal(1);
	std::array<int,2> nblocks({1,1});

	config::reader conf("config.json");
	int min_coarse = conf.get<int>("solver.min-coarse");

	do {
		ret.push_back(nblocks);
		if (nlx > nly) {
			if (((topoc.nglobal(0) / (nblocks[0]*2)) <= 2*min_coarse) or (npx / (nblocks[0]*2) <= 0))
				nblocks[1] *= 2;
			else
				nblocks[0] *= 2;
		} else {
			if (((topoc.nglobal(1) / (nblocks[1]*2)) <= 2*min_coarse) or (npy / (nblocks[1]*2) <=0))
				nblocks[0] *=2;
			else
				nblocks[1] *= 2;
		}

		nlx = topoc.nglobal(0) / nblocks[0];
		nly = topoc.nglobal(1) / nblocks[1];
	} while (keep_refining(npx, npy, nblocks[0], nblocks[1], nlx, nly, min_coarse));

	return ret;
}


namespace boxmg {
	bool operator<(const perf_state & s1, const perf_state & s2)
	{
		// TODO: handle this recursively
		// if (s1.model == nullptr) return true;
		// else if (s2.model == nullptr) return false;

		// len_t ng[2][2];
		// len_t nl[2][2];

		// for (auto i = 0; i < 2; i++) {
		// 	ng[0][i] = s1.model->grid(-1).nglobal(i);
		// 	ng[1][i] = s2.model->grid(-1).nglobal(i);
		// 	nl[0][i] = s1.model->grid(-1).nlocal(i);
		// 	nl[1][i] = s2.model->grid(-1).nlocal(i);
		// }

		// if (ng[0][0] < ng[1][0]) return true;
		// else if (ng[0][0] == ng[1][0] and ng[0][1] < ng[1][1]) return true;
		// else if (ng[0][0] == ng[1][0] and ng[0][1] == ng[1][1] and nl[0][0] < nl[1][0]) return true;
		// else if (ng[0][0] == ng[1][0] and ng[0][1] == ng[1][1] and nl[0][0] == nl[1][0] and nl[0][1] < nl[1][1]) return true;
		// else if
		return true;

	}
}
