#ifndef CEDAR_TIMELOG_H
#define CEDAR_TIMELOG_H

#include <limits>
#include <functional>
#include <string>
#include <vector>
#include <map>
#include <type_traits>
#include <chrono>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <cedar/util/log.h>
#include <cedar/mpi/redist_comms.h>
#include <cedar/util/timer_util.h>

namespace cedar {

	template<machine_mode mmode>
	class time_log
	{
	public:
		using timing_t = typename std::conditional<mmode==machine_mode::SERIAL,
			std::chrono::high_resolution_clock::time_point, double>::type;
		time_log(): lvl(0)
		{
			stimes.resize(1);
			ltimes.resize(1);
			counts.resize(1);
		}


		void begin(std::string label) {
			stimes[lvl][label] = timer_util<mmode>::now();
		}


		void end(std::string label) {
			auto endtime = timer_util<mmode>::now();
			auto elapsed = timer_util<mmode>::duration(stimes[lvl][label], endtime);
			ltimes[lvl][label] += elapsed;
			counts[lvl][label] += 1;
		}
		void up() { lvl--; }
		void down() {
			lvl++;
			if (static_cast<unsigned int>(lvl)+1 > stimes.size()) {
				stimes.resize(stimes.size()+1);
				ltimes.resize(ltimes.size()+1);
				counts.resize(counts.size()+1);
			}
		}

		void init(MPI_Comm comm) { this->comm = comm; }
		void redist(redist_comms comm) {rcomms.push_back(comm);redist_levels.push_back(lvl);}
		void save(std::string fname);

	private:
		MPI_Comm comm;
		std::vector<int> redist_levels;
		std::vector<redist_comms> rcomms;
		std::vector<std::map<std::string, double>> ltimes;     /** local elapsed times */
		std::vector<std::map<std::string, timing_t>> stimes;   /** start times */
		std::vector<std::map<std::string, int>> counts;
		int lvl;
	};
	template<> void time_log<machine_mode::MPI>::save(std::string fname);
	template<> void time_log<machine_mode::SERIAL>::save(std::string fname);

	extern time_log<machine_mode::MPI> tlog;
	extern time_log<machine_mode::SERIAL> tlog_ser;
	extern bool serial_timers;
	extern bool active;

	inline void timer_begin(std::string label) {
		if (active) {
			if (serial_timers)
				tlog_ser.begin(label);
			else
				tlog.begin(label);
		}
	}
	inline void timer_end(std::string label) {
		if (active) {
			if (serial_timers)
				tlog_ser.end(label);
			else
				tlog.end(label);
		}
	}
	inline void timer_up() {
		if (active) {
			if (serial_timers)
				tlog_ser.up();
			else
				tlog.up();
		}
	}
	inline void timer_down(){
		if (active) {
			if (serial_timers)
				tlog_ser.down();
			else
				tlog.down();
		}
	}
	inline void timer_pause() { active = false; }
	inline void timer_play() { active = true; }
	inline void timer_redist(redist_comms comm) { tlog.redist(comm); }
	inline void timer_init(MPI_Comm comm){ tlog.init(comm); serial_timers = false; }
	inline void timer_save(std::string fname){
		if (serial_timers)
			tlog_ser.save(fname);
		else
			tlog.save(fname);
	}
}

#endif
