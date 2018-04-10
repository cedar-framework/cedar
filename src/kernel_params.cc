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

	return params;
}


void ml_relax_params::init(config::reader & conf)
{
	this->enabled = conf.get<bool>("solver.ml-relax.enabled", false);
	this->min_gsz = conf.get<int>("solver.ml-relax.min-gsz", 3);
}

}
