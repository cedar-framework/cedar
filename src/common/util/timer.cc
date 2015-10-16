#include <mpi.h>

#include "log.h"
#include "timer.h"

using namespace boxmg;

void Timer::begin()
{
	if (log::timer.active())
		starttime = MPI_Wtime();
}


void Timer::end()
{
	if (log::timer.active()) {
		endtime = MPI_Wtime();
		log::timer << name << " time: " << (endtime - starttime) << std::endl;
	}
}
