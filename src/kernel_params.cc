#include <cedar/kernel_params.h>
#include <iostream>

using namespace cedar;

namespace cedar
{
std::shared_ptr<kernel_params> build_kernel_params(config::reader & conf)
{

	auto params = std::make_shared<kernel_params>();

	for (std::size_t i = 0; i < 3; ++i) params->periodic[i] = false;
	auto pers = conf.get<std::vector<bool>>("grid.periodic", {false,false,false});
	for (std::size_t i = 0; i < pers.size(); ++i) params->periodic[i] = pers[i];
	params->relax_symmetric = true;
	params->definite = true;
	params->ml_relax.init(conf);
	params->plane_config = conf.getconf("plane-config");
	if (params->plane_config == nullptr) {
		params->plane_config = std::make_shared<config::reader>("");
		params->plane_config->set("solver.relaxation", "line-xy");
		params->plane_config->set("solver.max-iter", 1);
	}

	return params;
}


void ml_relax_params::init(config::reader & conf)
{
	this->enabled = conf.get<bool>("solver.ml-relax.enabled", false);
	this->min_gsz = conf.get<int>("solver.ml-relax.min-gsz", 3);
}

}
