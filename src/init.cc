#include <cedar/global_manager.h>
#include <cedar/types.h>
#include <cedar/log.h>

namespace cedar {

static void wraplog()
{
	log::status.wrap(&gman.get<log::logger>().status);
	log::error.wrap(&gman.get<log::logger>().error);
	log::info.wrap(&gman.get<log::logger>().info);
	log::debug.wrap(&gman.get<log::logger>().debug);
	log::timer.wrap(&gman.get<log::logger>().timer);
	log::memory.wrap(&gman.get<log::logger>().memory);
}


void init(config & conf, MPI_Comm comm)
{
	auto params = build_global_params(conf);
	gman = global_manager<reg_globals>(params);
	gman.get<gmant::logger>().status << params;
	wraplog();
	gman.get<gmant::timer>().init(comm);
}

void init(config & conf)
{
	auto params = build_global_params(conf);
	gman = global_manager<reg_globals>(params);
	gman.get<gmant::logger>().status << params;
	wraplog();
}

}
