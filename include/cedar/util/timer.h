#ifndef CEDAR_TIMER_H
#define CEDAR_TIMER_H

#include <mpi.h>
#include <string>
#include <cedar/util/time_log.h>

namespace cedar {

	class timer
	{
	public:
		template<class T>
			timer(T&& name): name(std::forward<T>(name)) {};
		virtual void begin();
		virtual void end();
		virtual double time();
	protected:
		std::string name;
		double starttime, endtime;
	};

	class sync_timer : public timer
	{
	public:
		template<class T>
			sync_timer(MPI_Comm comm, T&& name) : timer(std::forward<T>(name)), comm(comm) {}
		virtual void begin();
		virtual void end();
	protected:
		MPI_Comm comm;
	};
}

#endif
