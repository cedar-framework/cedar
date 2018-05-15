#include <limits>
#include <random>

#include <cedar/perf/redist_generator.h>
#include <cedar/perf/greedy_iterator.h>
#include <cedar/perf/full_iterator.h>
#include <cedar/2d/util/topo.h>
#include <cedar/3d/util/topo.h>
#include <cedar/config.h>
#include <cedar/perf/cholesky_model.h>
#include <cedar/perf/const_model.h>
#include <cedar/perf/util.h>

#include <cedar/ss/astar.h>
#include <cedar/perf/search.h>

#include <cedar/perf/perf_factory.h>


using namespace cedar;

std::shared_ptr<vcycle_model> perf_factory::produce_vcycle(redist_settings & settings, int npx, int npy, len_t nx, len_t ny)
{
	using namespace cedar::cdr2;

	auto np = npx*npy;

	int min_coarse = settings.min_coarse;
	float tw = (1.0/(settings.machine_bandwidth)) * sizeof(real_t);
	float ts = settings.machine_latency;
	float tc = settings.machine_fprate;
	auto model = std::make_shared<vcycle_model>(2);
	model->set_comp_param(tc);
	model->set_comm_param(ts, tw);
	auto topo = util::model_topo(npx, npy, nx, ny);

	// auto nlx = topo->nlocal(0);
	// auto nly = topo->nlocal(1);

	int nlevels = compute_nlevels<2>(*topo, min_coarse);

	for (auto i = 0; i < nlevels; i++) {
		model->add_level(topo);
		topo = util::coarsen_topo(model->grid_ptr(0));
		// nlx = topo->nlocal(0);
		// nly = topo->nlocal(1);
	}

	if (np == 1) {
		auto & topoc = model->grid(0);
		auto cg_model = std::make_shared<cholesky_model>(topoc.nglobal(0)*topoc.nglobal(1));
		cg_model->set_comp_param(tc);
		model->set_cgperf(cg_model);
	}

	return model;
}


std::vector<std::vector<int>> get_choices(redist_settings & settings, int npx, int npy, len_t nx, len_t ny)
{
	using namespace cedar::cdr2;

	auto np = npx*npy;

	int min_coarse = settings.min_coarse;
	auto model = std::make_shared<vcycle_model>(2);

	auto topo = util::model_topo(npx, npy, nx, ny);

	// auto nlx = topo->nlocal(0);
	// auto nly = topo->nlocal(1);

	int nlevels = compute_nlevels<2>(*topo, min_coarse);

	for (auto i = 0; i < nlevels; i++) {
		model->add_level(topo);
		topo = util::coarsen_topo(model->grid_ptr(0));
		// nlx = topo->nlocal(0);
		// nly = topo->nlocal(1);
	}

	auto & topoc = model->grid(0);
	if (np == 1) {
		return {{0}};
	} else {
		// predict the best number of processor blocks
		std::vector<std::vector<int>> choices;
		int choice_num = 0;
		for (auto nblocks : redist_generator<full_iterator>({npx, npy}, {topoc.nglobal(0), topoc.nglobal(1)}, min_coarse)) {
			for (auto i : range(2))
				model->nblocks(i) = nblocks[i];
			auto paths = get_choices(settings, model->nblocks(0), model->nblocks(1),
			                         topoc.nglobal(0), topoc.nglobal(1));
			for (auto path : paths) {
				path.push_back(choice_num);
				choices.push_back(path);
			}
			choice_num++;
		}
		return choices;
	}
}


std::shared_ptr<vcycle_model> perf_factory::random_vcycle(redist_settings & settings, int npx, int npy, len_t nx, len_t ny, std::vector<int> path)
{
	using namespace cedar::cdr2;

	if (path.size() == 0) {
		auto choices = get_choices(settings, npx, npy, nx, ny);
		std::random_device rd;
		std::mt19937 rng(rd());
		std::uniform_int_distribution<int> uni(0, choices.size()-1);
		auto rand_choice = uni(rng);
		path = choices[rand_choice];
	}

	auto np = npx*npy;

	int min_coarse = settings.min_coarse;
	float tw = (1.0/(settings.machine_bandwidth)) * sizeof(real_t);
	float ts = settings.machine_latency;
	float tc = settings.machine_fprate;
	auto model = std::make_shared<vcycle_model>(2);
	model->set_comp_param(tc);
	model->set_comm_param(ts, tw);
	auto topo = util::model_topo(npx, npy, nx, ny);

	// auto nlx = topo->nlocal(0);
	// auto nly = topo->nlocal(1);

	int nlevels = compute_nlevels<2>(*topo, min_coarse);

	for (auto i = 0; i < nlevels; i++) {
		model->add_level(topo);
		topo = util::coarsen_topo(model->grid_ptr(0));
		// nlx = topo->nlocal(0);
		// nly = topo->nlocal(1);
	}

	auto & topoc = model->grid(0);
	if (np == 1) {
		auto cg_model = std::make_shared<cholesky_model>(topoc.nglobal(0)*topoc.nglobal(1));
		cg_model->set_comp_param(tc);
		model->set_cgperf(cg_model);
	} else {
		// predict the best number of processor blocks
		std::vector<std::array<int,2>> choices;
		auto redist_subsets = redist_generator<full_iterator>({npx, npy}, {topoc.nglobal(0), topoc.nglobal(1)}, min_coarse);
		for (auto nblocks : redist_subsets) {
			for (auto i : range(2))
				model->nblocks(i) = nblocks[i];

			choices.push_back(std::array<int,2>({model->nblocks(0), model->nblocks(1)}));
		}

		auto rand_choice = path.back();
		path.pop_back();

		model->nblocks(0) = choices[rand_choice][0];
		model->nblocks(1) = choices[rand_choice][1];

		auto cg_model = perf_factory::random_vcycle(settings, model->nblocks(0), model->nblocks(1),
		                                            topoc.nglobal(0), topoc.nglobal(1), path);
		model->set_cgperf(cg_model);
	}

	return model;
}


std::shared_ptr<vcycle_model> perf_factory::manual_vcycle(redist_settings & settings, int npx, int npy, len_t nx, len_t ny)
{
	using namespace cedar::cdr2;

	auto np = npx*npy;

	auto model = perf_factory::produce_vcycle(settings, npx, npy, nx, ny);

	if (np > 1) {
		auto & path = settings.path;
		bool found = false;
		for (auto & proc : path) {
			if (found) {
				model->nblocks(0) = proc[0];
				model->nblocks(1) = proc[1];
				break;
			} else if (proc[0] == npx and proc[1] == npy) {
				found = true;
			}
		}

		auto cg_model = perf_factory::manual_vcycle(settings, model->nblocks(0), model->nblocks(1),
		                                            model->grid(0).nglobal(0), model->grid(0).nglobal(1));
		model->set_cgperf(cg_model);
	}

	return model;
}


std::shared_ptr<vcycle_model> perf_factory::dfs_vcycle(redist_settings & settings, int npx, int npy, len_t nx, len_t ny, bool terminate, int rlevel)
{
	using namespace cedar::cdr2;

	auto np = npx*npy;

	int min_coarse = settings.min_coarse;
	float tw = (1.0/(settings.machine_bandwidth)) * sizeof(real_t);
	float ts = settings.machine_latency;
	float tc = settings.machine_fprate;
	auto model = std::make_shared<vcycle_model>(2);
	model->set_comp_param(tc);
	model->set_comm_param(ts, tw);
	auto topo = util::model_topo(npx, npy, nx, ny);

	// auto nlx = topo->nlocal(0);
	// auto nly = topo->nlocal(1);

	int nlevels = compute_nlevels<2>(*topo, min_coarse);

	for (auto i = 0; i < nlevels; i++) {
		model->add_level(topo);
		topo = util::coarsen_topo(model->grid_ptr(0));
		// nlx = topo->nlocal(0);
		// nly = topo->nlocal(1);
	}

	auto & topoc = model->grid(0);
	if (terminate and np != 1) {
		auto cg_model = perf_factory::dfs_vcycle(settings, 1,1, topoc.nglobal(0), topoc.nglobal(1), true);
		cg_model->set_comm_param(0, 0); // Since this is serial
		model->nblocks(0) = 1;
		model->nblocks(1) = 1;
		model->set_cgperf(cg_model);
	} else if (np == 1) {
		auto cg_model = std::make_shared<cholesky_model>(topoc.nglobal(0)*topoc.nglobal(1));
		cg_model->set_comp_param(tc);
		model->set_cgperf(cg_model);
	} else {
		// predict the best number of processor blocks
		std::array<int,2> best_blocks;
		float best_time = std::numeric_limits<float>::max();
		std::shared_ptr<vcycle_model> best_cg;
		for (auto nblocks : redist_generator<greedy_iterator>({npx, npy}, {topoc.nglobal(0), topoc.nglobal(1)},
		                                                      min_coarse))
		{
			model->nblocks(0) = nblocks[0];
			model->nblocks(1) = nblocks[1];
			auto cg_model = perf_factory::dfs_vcycle(settings, model->nblocks(0), model->nblocks(1),
			                                         topoc.nglobal(0), topoc.nglobal(1),false,rlevel+1);
			// set coarse solve to 0
			// cg_model->set_cgperf(std::make_shared<const_model>(0));

			model->set_cgperf(cg_model);
			float this_time = model->tcgsolve();
			if (this_time < best_time) {
				best_time = this_time;
				best_blocks[0] = model->nblocks(0);
				best_blocks[1] = model->nblocks(1);
				best_cg = cg_model;
			}
		}

		model->nblocks(0) = best_blocks[0];
		model->nblocks(1) = best_blocks[1];
		model->set_cgperf(best_cg);
	}

	return model;
}


std::shared_ptr<vcycle_model> perf_factory::astar_vcycle(redist_settings & settings, int npx, int npy, len_t nx, len_t ny)
{
	using namespace cedar;
	perf_problem pprob(settings);
	pprob.initial_state = perf_state();

	pprob.initial_state.model = perf_factory::produce_vcycle(settings,npx,npy,nx,ny);

	using node_ptr = std::shared_ptr<perf_node>;
	std::function<float(node_ptr)> heuristic = [](node_ptr nd) {
		auto model = nd->state.model;
		auto nlevels = std::log2(model->grid(0).nproc());
		return nlevels * model->tsmooth(-1);
	};
	auto sol = ss::astar<perf_solution>(pprob, heuristic);

	return sol.model();
}


std::array<len_t,2> perf_factory::graph_vcycle(std::ostream & os, int npx, int npy, len_t nx, len_t ny, bool terminate, int rlevel)
{
	// TODO: replace with graph library (libcgraph)
	using namespace cedar::cdr2;

	const bool FULL = false;

	if (rlevel == 0) {
		os << "digraph {" << '\n';
		os << "a [label=\"proc grid\ncoarse grid\"];" << '\n';
	}

	auto node_id = [](int rlevel, int npx, int npy, len_t ngx, len_t ngy, len_t ngcx, len_t ngcy) {
		if (FULL)
			return "a" + std::to_string(rlevel) + "_" + std::to_string(npx) + "_" + std::to_string(npy) + "_" +
				std::to_string(ngx) + "_" + std::to_string(ngy);
		else
			return "a_" + std::to_string(npx) + "_" + std::to_string(npy) + "_" +
				std::to_string(ngcx) + "_" + std::to_string(ngcy);
	};

	auto np = npx*npy;

	config conf("config.json");

	int min_coarse = conf.get<int>("solver.min-coarse");
	float tw = (1.0/conf.get<float>("machine.bandwidth")) * sizeof(real_t);
	float ts = conf.get<float>("machine.latency");
	float tc = conf.get<float>("machine.fp_perf");
	auto model = std::make_shared<vcycle_model>(2);
	model->set_comp_param(tc);
	model->set_comm_param(ts, tw);
	auto topo = util::model_topo(npx, npy, nx, ny);

	// auto nlx = topo->nlocal(0);
	// auto nly = topo->nlocal(1);

	int nlevels = compute_nlevels<2>(*topo, min_coarse);

	for (auto i = 0; i < nlevels; i++) {
		model->add_level(topo);
		topo = util::coarsen_topo(model->grid_ptr(0));
		// nlx = topo->nlocal(0);
		// nly = topo->nlocal(1);
	}

	auto & topoc = model->grid(0);

	os << node_id(rlevel, npx, npy, nx, ny, topoc.nglobal(0), topoc.nglobal(1)) << " " << "[label=\""
	   << npx << " x " << npy << "\\n"
	   << topoc.nglobal(0)-2 << " x " << topoc.nglobal(1)-2 << "\"];\n";

	if (terminate and np != 1) {
		auto ngc = perf_factory::graph_vcycle(os, 1,1, topoc.nglobal(0), topoc.nglobal(1), true, rlevel+1);
		os << node_id(rlevel, npx, npy, nx, ny, topoc.nglobal(0), topoc.nglobal(1)) << " -> "
		   << node_id(rlevel+1, model->nblocks(0), model->nblocks(1), topoc.nglobal(0), topoc.nglobal(1), ngc[0], ngc[1]) << ";\n";
	} else if (np == 1) {
		auto cg_model = std::make_shared<cholesky_model>(topoc.nglobal(0)*topoc.nglobal(1));
		cg_model->set_comp_param(tc);
		model->set_cgperf(cg_model);
	} else {
		// predict the best number of processor blocks
		for (auto nblocks : redist_generator<greedy_iterator>({npx, npy}, {topoc.nglobal(0), topoc.nglobal(1)}, min_coarse)) {
			for (auto i : range(2))
				model->nblocks(i) = nblocks[i];

			auto ngc = perf_factory::graph_vcycle(os, model->nblocks(0), model->nblocks(1),
			                                      topoc.nglobal(0), topoc.nglobal(1),false,rlevel+1);
			os << node_id(rlevel, npx, npy, nx, ny, topoc.nglobal(0), topoc.nglobal(1)) << " -> "
			   << node_id(rlevel+1, model->nblocks(0), model->nblocks(1), topoc.nglobal(0), topoc.nglobal(1), ngc[0], ngc[1]) << ";\n";
		}
	}

	if (rlevel == 0) {
		os << "}";
	}

	std::array<len_t, 2> ret({topoc.nglobal(0), topoc.nglobal(1)});

	return ret;
}

// std::shared_ptr<vcycle_model> perf_factory::produce_vcycle(int npx, int npy, int npz,
//                                                            len_t nx, len_t ny, len_t nz, bool terminate)
// {
// 	using namespace cedar::cdr3;

// 	config conf("perf.json");
// 	int min_coarse = conf.get<int>("solver.min-coarse", 3);
// 	auto model = std::make_shared<vcycle_model>(3);
// 	auto topo = util::model_topo(npx, npy, npz, nx, ny, nz);

// 	auto nlx = topo->nlocal(0);
// 	auto nly = topo->nlocal(1);
// 	auto nlz = topo->nlocal(2);

// 	do {
// 		model->add_level(topo);
// 		topo = util::coarsen_topo(model->grid_ptr(0));
// 		nlx = topo->nlocal(0);
// 		nly = topo->nlocal(1);
// 		nlz = topo->nlocal(2);
// 	} while (std::min({nlx, nly, nlz}) >= min_coarse);

// 	auto & topoc = model->grid(0);
// 	if (np != 1) {
// 		auto cg_model = perf_factory::produce_vcycle(1, topoc.nglobal(0), topoc.nglobal(1), topoc.nglobal(2));
// 		model->set_cgperf(cg_model);
// 	} else {
// 		auto cg_model = std::make_shared<cholesky_model>(topoc.nglobal(0)*topoc.nglobal(1)*topoc.nglobal(2));
// 		model->set_cgperf(cg_model);
// 	}

// 	return model;
// }
