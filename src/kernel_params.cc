#include <cedar/kernel_params.h>
#include <cedar/util/log.h>
#include <iostream>

using namespace cedar;

namespace cedar
{

std::ostream & operator<<(std::ostream & os, const kernel_params & obj)
{
	os << '\n';
	os << "---------------\n";
	os << "Kernel Settings\n";
	os << "---------------\n";
	os << "symmetric relaxation: ";
	if (obj.relax_symmetric)
		os << "yes\n";
	else
		os << "no\n";

	os << "periodic x:           ";
	if (obj.periodic[0])
		os << "yes\n";
	else
		os << "no\n";
	os << "periodic y:           ";
	if (obj.periodic[1])
		os << "yes\n";
	else
		os << "no\n";
	os << "periodic z:           ";
	if (obj.periodic[2])
		os << "yes\n";
	else
		os << "no\n";

	os << "multilevel lines:     ";
	if (obj.ml_relax.enabled) {
		os << "yes\n";
		os << "ml-relax group size:  ";
		os << obj.ml_relax.min_gsz << '\n';
		os << "ml-relax factorize :  ";
		if (obj.ml_relax.factorize)
			os << "yes\n";
		else
			os << "no\n";
	} else
		os << "no\n";

	os << "halo exchange:        ";
	if (obj.halo == kernel_params::halo_lib::msg)
		os << "msg";
	else
		os << "tausch";

	return os;
}


std::shared_ptr<kernel_params> build_kernel_params(config & conf)
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
		params->plane_config = std::make_shared<config>("");
		params->plane_config->set("solver.relaxation", "line-xy");
		params->plane_config->set("solver.max-iter", 1);
		params->plane_config->set("halo-exchange", "tausch");
	}

	auto halo_exchange_name = conf.get<std::string>("halo-exchange", "msg");
	if (halo_exchange_name == "msg")
		params->halo = kernel_params::halo_lib::msg;
	else if (halo_exchange_name == "tausch")
		params->halo = kernel_params::halo_lib::tausch;
	else
		log::error << "invalid halo-exchange: " << halo_exchange_name << std::endl;

	return params;
}


void ml_relax_params::init(config & conf)
{
	this->enabled = conf.get<bool>("solver.ml-relax.enabled", false);
	this->min_gsz = conf.get<int>("solver.ml-relax.min-gsz", 3);
	this->factorize = conf.get<bool>("solver.ml-relax.factorize", true);
}

}
