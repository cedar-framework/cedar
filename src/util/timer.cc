#include <cedar/log.h>
#include <cedar/util/timer.h>

using namespace cedar;


void timer::begin()
{
	if (log::timer.active())
		starttime = MPI_Wtime();
}


void timer::end()
{
	if (log::timer.active()) {
		endtime = MPI_Wtime();
		log::timer << name << " time: " << (endtime - starttime) << std::endl;
	}
}


double timer::time()
{
	return endtime - starttime;
}


void sync_timer::begin()
{
	if (log::timer.active())
		MPI_Barrier(comm);
	timer::begin();
}


void sync_timer::end()
{
	if (log::timer.active())
		MPI_Barrier(comm);
	timer::end();
}

