#include <cedar/util/log.h>
#include <cedar/capi.h>
#include <cedar/interface/object.h>
#include <cedar/interface/config.h>

extern "C"
{
	int cedar_log_init(cedar_config conf)
	{
		auto confobj = cedar_config_getobj(conf);
		if (not confobj)
			return CEDAR_ERR_CONFIG;

		cedar::log::init(*confobj);

		return CEDAR_SUCCESS;
	}
}
