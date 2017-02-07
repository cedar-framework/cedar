#include "boxmg/util/log.h"
#include <boxmg/util/time_log.h>

extern "C" {

	void print_error(char *string)
	{
		boxmg::log::error << string << std::endl;
	}

	void ftimer_begin(char *string)
	{
		boxmg::timer_begin(string);
	}

	void ftimer_end(char *string)
	{
		boxmg::timer_end(string);
	}
}
