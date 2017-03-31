#include <cedar/interface/c/timer.h>
#include <cedar/util/time_log.h>

extern "C"
{

	void bmg_timer_save(const char * fname)
	{
		std::string fname_cpp(fname);
		cedar::timer_save(fname_cpp);
	}
}
