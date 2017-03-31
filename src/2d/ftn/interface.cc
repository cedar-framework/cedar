#include "cedar/util/log.h"
#include <cedar/util/time_log.h>

extern "C" {

	void print_error(char *string)
	{
		cedar::log::error << string << std::endl;
	}

	void ftimer_begin(char *string)
	{
		cedar::timer_begin(string);
	}

	void ftimer_end(char *string)
	{
		cedar::timer_end(string);
	}
}
