#include "cedar/util/log.h"
#include <cedar/types.h>
#include <cedar/util/time_log.h>
#include <cedar/halo_exchanger.h>

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

	void halo_exchange(int k, int nog, cedar::real_t *q, cedar::len_t II,
	                   cedar::len_t JJ, cedar::halo_exchanger<2> *halof, void *ctx)
	{
		halof->exchange(k, nog, q, {II, JJ}, ctx);
	}
}
