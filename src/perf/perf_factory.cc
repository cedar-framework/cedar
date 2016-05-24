#include <limits>
#include <boxmg/2d/util/topo.h>
#include <boxmg/3d/util/topo.h>
#include <boxmg/config/reader.h>
#include <boxmg/perf/cholesky_model.h>
#include <boxmg/perf/const_model.h>
#include <boxmg/perf/params.h>
#include <boxmg/perf/util.h>

#include <boxmg/ss/astar.h>
#include <boxmg/perf/search.h>

#include <boxmg/perf/perf_factory.h>


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


std::shared_ptr<vcycle_model> perf_factory::produce_vcycle(int npx, int npy, len_t nx, len_t ny)
{
	using namespace boxmg::bmg2d;

	auto np = npx*npy;

	config::reader conf("config.json");

	int min_coarse = conf.get<int>("solver.min-coarse");
	float ts = conf.get<float>("machine.bandwidth");
	float tw = conf.get<float>("machine.latency");
	float tc = conf.get<float>("machine.fp_perf");
	auto model = std::make_shared<vcycle_model>(2);
	model->set_comp_param(params::compute_tc(2, conf));
	model->set_comm_param(ts, tw);
	auto topo = util::model_topo(npx, npy, nx, ny);

	auto nlx = topo->nlocal(0);
	auto nly = topo->nlocal(1);

	int nlevels = compute_nlevels<2>(*topo, min_coarse);

	for (auto i = 0; i < nlevels; i++) {
		model->add_level(topo);
		topo = util::coarsen_topo(model->grid_ptr(0));
		nlx = topo->nlocal(0);
		nly = topo->nlocal(1);
	}

	if (np == 1) {
		auto & topoc = model->grid(0);
		auto cg_model = std::make_shared<cholesky_model>(topoc.nglobal(0)*topoc.nglobal(1));
		cg_model->set_comp_param(tc);
		model->set_cgperf(cg_model);
	}

	return model;
}


std::shared_ptr<vcycle_model> perf_factory::dfs_vcycle(int npx, int npy, len_t nx, len_t ny, bool terminate, int rlevel)
{
	using namespace boxmg::bmg2d;

	auto np = npx*npy;

	config::reader conf("config.json");

	int min_coarse = conf.get<int>("solver.min-coarse");
	float ts = conf.get<float>("machine.bandwidth");
	float tw = conf.get<float>("machine.latency");
	float tc = conf.get<float>("machine.fp_perf");
	auto model = std::make_shared<vcycle_model>(2);
	model->set_comp_param(params::compute_tc(2, conf));
	model->set_comm_param(ts, tw);
	auto topo = util::model_topo(npx, npy, nx, ny);

	auto nlx = topo->nlocal(0);
	auto nly = topo->nlocal(1);

	int nlevels = compute_nlevels<2>(*topo, min_coarse);

	for (auto i = 0; i < nlevels; i++) {
		model->add_level(topo);
		topo = util::coarsen_topo(model->grid_ptr(0));
		nlx = topo->nlocal(0);
		nly = topo->nlocal(1);
	}

	auto & topoc = model->grid(0);
	if (terminate and np != 1) {
		auto cg_model = perf_factory::dfs_vcycle(1,1, topoc.nglobal(0), topoc.nglobal(1), true);
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
		model->nblocks(0) = 1;
		model->nblocks(1) = 1;
		std::array<int,2> best_blocks;
		float best_time = std::numeric_limits<float>::max();
		std::shared_ptr<vcycle_model> best_cg;
		len_t nlx = topoc.nglobal(0);
		len_t nly = topoc.nglobal(1);
		do {
			auto cg_model = perf_factory::dfs_vcycle(model->nblocks(0), model->nblocks(1),
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
			//if ((npx / model->nblocks(0)) > (npy / model->nblocks(1))) {
			// greedily adding processor blocks by local problem size
			if (nlx > nly) {
				if (((topoc.nglobal(0) / (model->nblocks(0)*2)) <= 2*min_coarse) or (npx / (model->nblocks(0)*2) <= 0))
					model->nblocks(1) *= 2;
				else
					model->nblocks(0) *= 2;
			} else {
				if (((topoc.nglobal(1) / (model->nblocks(1)*2)) <= 2*min_coarse) or (npy / (model->nblocks(1)*2) <=0))
					model->nblocks(0) *=2;
				else
					model->nblocks(1) *= 2;
			}

			nlx = topoc.nglobal(0) / model->nblocks(0);
			nly = topoc.nglobal(1) / model->nblocks(1);
		} while (keep_refining(npx, npy, model->nblocks(0), model->nblocks(1),
		                       nlx, nly, min_coarse));


		model->nblocks(0) = best_blocks[0];
		model->nblocks(1) = best_blocks[1];
		model->set_cgperf(best_cg);
	}

	return model;
}


std::shared_ptr<vcycle_model> perf_factory::astar_vcycle(int npx, int npy, len_t nx, len_t ny)
{
	using namespace boxmg;
	perf_problem pprob;
	pprob.initial_state = perf_state();

	pprob.initial_state.model = perf_factory::produce_vcycle(npx,npy,nx,ny);

	using node_ptr = std::shared_ptr<perf_node>;
	auto heuristic = [](node_ptr nd) {
		auto model = nd->state.model;
		auto topo = bmg2d::util::model_topo(1, 1, model->grid(0).nglobal(0), model->grid(0).nglobal(1));
		int nlevels = compute_nlevels<2>(*topo, 3);

		return model->tsmooth(-1)*nlevels;
	};
	auto sol = ss::astar<perf_solution>(pprob, heuristic);

	return sol.model();
}

// std::shared_ptr<vcycle_model> perf_factory::produce_vcycle(int npx, int npy, int npz,
//                                                            len_t nx, len_t ny, len_t nz, bool terminate)
// {
// 	using namespace boxmg::bmg3;

// 	config::reader conf("perf.json");
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
