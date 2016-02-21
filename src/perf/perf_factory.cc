#include <boxmg/2d/util/topo.h>
#include <boxmg/3d/util/topo.h>
#include <boxmg/config/reader.h>
#include <boxmg/perf/cholesky_model.h>
#include <boxmg/perf/params.h>

#include <boxmg/perf/perf_factory.h>


using namespace boxmg;

std::shared_ptr<vcycle_model> perf_factory::produce_vcycle(int np, len_t nx, len_t ny)
{
	using namespace boxmg::bmg2d;

	config::Reader conf("perf.json");

	int min_coarse = conf.get<int>("solver.min-coarse", 3);
	float ts = conf.get<float>("machine.bandwidth", 1e-6);
	float tw = conf.get<float>("machine.latency", 1e-4);
	float tc = conf.get<float>("machine.fp_perf", 1e-8);
	auto model = std::make_shared<vcycle_model>(2);
	model->set_comp_param(params::compute_tc(2));
	model->set_comm_param(ts, tw);
	auto topo = util::create_topo(np, nx, ny);

	auto nlx = topo->nlocal(0);
	auto nly = topo->nlocal(1);

	do {
		model->add_level(topo);
		topo = util::coarsen_topo(model->grid_ptr(0));
		nlx = topo->nlocal(0);
		nly = topo->nlocal(1);
		log::status << nlx << " " << nly << std::endl;
	} while (std::min(nlx, nly) >= min_coarse);

	auto & topoc = model->grid(0);
	if (np != 1) {
		log::status << np << std::endl;
		auto cg_model = perf_factory::produce_vcycle(1, topoc.nglobal(0), topoc.nglobal(1));
		cg_model->set_comm_param(0, 0); // Since this is serial
		model->set_cgperf(cg_model);
	} else {
		auto cg_model = std::make_shared<cholesky_model>(topoc.nglobal(0)*topoc.nglobal(1));
		model->set_cgperf(cg_model);
		model->set_comp_param(tc);
	}

	return model;
}


std::shared_ptr<vcycle_model> perf_factory::produce_vcycle(int np, len_t nx, len_t ny, len_t nz)
{
	using namespace boxmg::bmg3;

	config::Reader conf("perf.json");
	int min_coarse = conf.get<int>("solver.min-coarse", 3);
	auto model = std::make_shared<vcycle_model>(3);
	auto topo = util::create_topo(np, nx, ny, nz);

	auto nlx = topo->nlocal(0);
	auto nly = topo->nlocal(1);
	auto nlz = topo->nlocal(2);

	do {
		model->add_level(topo);
		topo = util::coarsen_topo(model->grid_ptr(0));
		nlx = topo->nlocal(0);
		nly = topo->nlocal(1);
		nlz = topo->nlocal(2);
	} while (std::min({nlx, nly, nlz}) >= min_coarse);

	auto & topoc = model->grid(0);
	if (np != 1) {
		auto cg_model = perf_factory::produce_vcycle(1, topoc.nglobal(0), topoc.nglobal(1), topoc.nglobal(2));
		model->set_cgperf(cg_model);
	} else {
		auto cg_model = std::make_shared<cholesky_model>(topoc.nglobal(0)*topoc.nglobal(1)*topoc.nglobal(2));
		model->set_cgperf(cg_model);
	}

	return model;
}
