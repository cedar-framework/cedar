#ifndef CEDAR_PERF_SEARCH_H
#define CEDAR_PERF_SEARCH_H

#include <array>
#include <vector>

#include <cedar/ss/node.h>
#include <cedar/ss/problem.h>
#include <cedar/ss/solution.h>
#include <cedar/perf/vcycle_model.h>
#include <cedar/perf/cholesky_model.h>

namespace cedar {

struct perf_state
{
	std::shared_ptr<vcycle_model> model;
	friend bool operator<(const perf_state & s1, const perf_state & s2);
};


using perf_node = ss::node<perf_state, std::array<int, 2>, float>;

struct perf_problem : ss::problem<perf_state, std::array<int,2>, float>
{
    perf_problem(config::reader & conf): conf(conf) {}
	bool goal_test(perf_state & model);
	perf_state result(perf_state & model, std::array<int,2> action);
	float step_cost(perf_state & state, std::array<int,2> act);
	std::vector<std::array<int,2>> actions(perf_state & state);

private:
	config::reader & conf;
};


class perf_solution : public ss::solution<perf_state, std::array<int,2>, float>
{
public:
	using node_ptr = std::shared_ptr<perf_node>;
perf_solution(node_ptr sol_node) :
	ss::solution<perf_state, std::array<int,2>,float>(sol_node)
	{
	}

	std::shared_ptr<vcycle_model> model() {
		auto pnode = snode();
		auto model = pnode->state.model;
		auto & topoc = model->grid(0);
		auto cg_model = std::make_shared<cholesky_model>(topoc.nglobal(0)*topoc.nglobal(1));
		cg_model->set_comp_param(model->get_comp_param());
		model->set_cgperf(cg_model);

		while (pnode->parent != nullptr) {
			pnode->parent->state.model->set_cgperf(pnode->state.model);
			pnode->parent->state.model->nblocks(0) = pnode->state.model->grid(0).nproc(0);
			pnode->parent->state.model->nblocks(1) = pnode->state.model->grid(0).nproc(1);
			pnode = pnode->parent;
		}

		return pnode->state.model;
	}
};

}

#endif
