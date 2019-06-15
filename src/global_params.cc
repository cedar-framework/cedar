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

global_params build_global_params(config & conf)
{
	global_params params;

	params.log_level = getloglevel(conf);
	auto mtype = conf.get<std::string>("memory", "system");
	if (mtype == "managed")
		params.memory_type = memtype::managed;
	else
		params.memory_type = memtype::system;

	return params;
}


std::ostream & operator<<(std::ostream & os, const global_params & obj)
{
	os << '\n';
	os << "---------------\n";
	os << "Global Settings\n";
	os << "---------------\n";

	os << "memory type: ";
	if (obj.memory_type == memtype::managed)
		os << "managed\n";
	else
		os << "system\n";

	os << "log levels:  ";
	if (obj.log_level & loglevel::status)
		os << "status ";
	if (obj.log_level & loglevel::info)
		os << "info ";
	if (obj.log_level & loglevel::error)
		os << "error ";
	if (obj.log_level & loglevel::memory)
		os << "memory ";
	if (obj.log_level & loglevel::debug)
		os << "debug ";
	if (obj.log_level & loglevel::timer)
		os << "timer";
	os << "\n";

	return os;
}
}
