#include <map>

#include <cedar/global_params.h>

namespace cedar {

loglevel_t getloglevel(config & conf)
{
	std::map<std::string, loglevel_t> lmap {
		{"status", loglevel::status},
		{"info"  , loglevel::info},
		{"error" , loglevel::error},
		{"memory", loglevel::memory},
		{"debug" , loglevel::debug},
		{"timer",  loglevel::timer}};

	loglevel_t level = 0;
	auto clevels = conf.getvec<std::string>("log");
	for (auto &clvl: clevels) level |= lmap[clvl];

	return level;
}

std::shared_ptr<global_params> build_global_params(config & conf)
{
	auto params = std::make_shared<global_params>();

	params->log_level = getloglevel(conf);

	return params;
}

}
