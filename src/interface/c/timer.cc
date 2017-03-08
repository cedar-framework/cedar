#include <boxmg/interface/c/timer.h>
#include <boxmg/util/time_log.h>

extern "C"
{

	void bmg_timer_save(const char * fname)
	{
		std::string fname_cpp(fname);
		boxmg::timer_save(fname_cpp);
	}
}
