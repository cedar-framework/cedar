#include <cedar/global_manager.h>
#include <cedar/types.h>

namespace cedar {

void init(config & conf)
{
	auto params = build_global_params(conf);
	gman = global_manager<reg_globals>(*params);
	log::status.wrap(&gman.get<log::logger>().status);
	log::error.wrap(&gman.get<log::logger>().error);
	log::info.wrap(&gman.get<log::logger>().info);
	log::debug.wrap(&gman.get<log::logger>().debug);
	log::timer.wrap(&gman.get<log::logger>().timer);
	log::memory.wrap(&gman.get<log::logger>().memory);
}

}
