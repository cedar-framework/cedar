#include <limits>
#include <boxmg/2d/util/topo.h>
#include <boxmg/3d/util/topo.h>
#include <boxmg/config/reader.h>
#include <boxmg/perf/cholesky_model.h>
#include <boxmg/perf/params.h>

#include <boxmg/perf/perf_factory.h>


using namespace boxmg;

std::shared_ptr<vcycle_model> perf_factory::produce_vcycle(int npx, int npy, len_t nx, len_t ny, bool terminate)
{
	using namespace boxmg::bmg2d;

	auto np = npx*npy;

	config::reader conf("perf.json");

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

	do {
		model->add_level(topo);
		topo = util::coarsen_topo(model->grid_ptr(0));
		nlx = topo->nlocal(0);
		nly = topo->nlocal(1);
	} while (std::min(nlx, nly) >= min_coarse);

	auto & topoc = model->grid(0);
	if (terminate and np != 1) {
		auto cg_model = perf_factory::produce_vcycle(1,1, topoc.nglobal(0), topoc.nglobal(1), true);
		cg_model->set_comm_param(0, 0); // Since this is serial
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
		while ((npx/model->nblocks(0)) > 0 and (npy/model->nblocks(1) > 0)) {
			auto cg_model = perf_factory::produce_vcycle(npx / model->nblocks(0), npy / model->nblocks(1),
			                                             topoc.nglobal(0), topoc.nglobal(1), true);
			model->set_cgperf(cg_model);
			float this_time = model->tcgsolve();
			log::error << "cgsolve with " << model->nblocks(0)*model->nblocks(1) << " chunks would take " << this_time << " seconds." << std::endl;
			if (this_time < best_time) {
				best_time = this_time;
				best_blocks[0] = model->nblocks(0);
				best_blocks[1] = model->nblocks(1);
			}
			if ((npx / model->nblocks(0)) > (npy / model->nblocks(1)))
				model->nblocks(0) *= 2;
			else
				model->nblocks(1) *= 2;
		}

		model->nblocks(0) = best_blocks[0];
		model->nblocks(1) = best_blocks[1];
		auto cg_model = perf_factory::produce_vcycle(npx / model->nblocks(0), npy / model->nblocks(1),
		                                             topoc.nglobal(0), topoc.nglobal(1));
		model->set_cgperf(cg_model);
	}

	return model;
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
