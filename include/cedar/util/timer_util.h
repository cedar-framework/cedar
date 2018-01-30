#ifndef CEDAR_TIMERUTIL_H
#define CEDAR_TIMERUTIL_H

#include <mpi.h>
#include <type_traits>
#include <chrono>

namespace cedar {

	enum class machine_mode {SERIAL, MPI};

	template<machine_mode mmode>
	struct timer_util
	{
		using timing_t = typename std::conditional<mmode==machine_mode::SERIAL,
		                                           std::chrono::high_resolution_clock::time_point,
		                                           double>::type;
		static double duration(timing_t beg, timing_t end);
		static timing_t now();
	};

	template<> inline double timer_util<machine_mode::SERIAL>::duration(timer_util<machine_mode::SERIAL>::timing_t beg,
	                                                           timer_util<machine_mode::SERIAL>::timing_t end) {
		return std::chrono::duration_cast<std::chrono::duration<double>>(end-beg).count();
	}

	template<> inline timer_util<machine_mode::SERIAL>::timing_t timer_util<machine_mode::SERIAL>::now() {
		return std::chrono::high_resolution_clock::now();
	}

	template<> inline double timer_util<machine_mode::MPI>::duration(timer_util<machine_mode::MPI>::timing_t beg,
	                                                        timer_util<machine_mode::MPI>::timing_t end) {
		return end - beg;
	}

	template<> inline timer_util<machine_mode::MPI>::timing_t timer_util<machine_mode::MPI>::now() { return MPI_Wtime(); }


}

#endif
