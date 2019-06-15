#ifndef CEDAR_UTIL_AGG_TIMER_H
#define CEDAR_UTIL_AGG_TIMER_H

#include <limits>
#include <functional>
#include <string>
#include <vector>
#include <map>
#include <type_traits>
#include <chrono>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <cedar/mpi/redist_comms.h>
#include <cedar/util/timer_util.h>
#include <cedar/global_params.h>

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
		active.resize(1);
	}


	void begin(std::string label) {
		if (not active[lvl][label]) {
			stimes[lvl][label] = timer_util<mmode>::now();
			active[lvl][label] = true;
		}
	}


	void end(std::string label) {
		if (active[lvl][label]) {
			auto endtime = timer_util<mmode>::now();
			auto elapsed = timer_util<mmode>::duration(stimes[lvl][label], endtime);
			ltimes[lvl][label] += elapsed;
			counts[lvl][label] += 1;
			active[lvl][label] = false;
		}
	}
	void up() { lvl--; }
	void down() {
		lvl++;
		if (static_cast<unsigned int>(lvl)+1 > stimes.size()) {
			stimes.resize(stimes.size()+1);
			ltimes.resize(ltimes.size()+1);
			counts.resize(counts.size()+1);
			active.resize(active.size()+1);
		}
	}

	void init(MPI_Comm comm) { this->comm = comm; }
	void redist(redist_comms comm) {rcomms.push_back(comm);redist_levels.push_back(lvl);}
	void save(const std::string & fname);

private:
	MPI_Comm comm;
	std::vector<int> redist_levels;
	std::vector<redist_comms> rcomms;
	std::vector<std::map<std::string, double>> ltimes;     /** local elapsed times */
	std::vector<std::map<std::string, timing_t>> stimes;   /** start times */
	std::vector<std::map<std::string, int>> counts;
	std::vector<std::map<std::string, bool>> active;
	int lvl;
};
template<> void time_log<machine_mode::MPI>::save(const std::string & fname);
template<> void time_log<machine_mode::SERIAL>::save(const std::string & fname);

class agg_timer
{
public:
	agg_timer(global_params & params);

	void begin(const std::string & label)
	{
		if (active) {
			if (serial)
				slog.begin(label);
			else
				plog.begin(label);
		}
	}

	void end(const std::string & label)
	{
		if (active) {
			if (serial)
				slog.end(label);
			else
				plog.end(label);
		}
	}

	void up()
	{
		if (active) {
			if (serial)
				slog.up();
			else
				plog.up();
		}
	}

	void down()
	{
		if (active) {
			if (serial)
				slog.down();
			else
				plog.down();
		}
	}

	void pause() { active = false; }
	void play() { active = true; }
	void redist(redist_comms comm) { plog.redist(comm); }
	void init(MPI_Comm comm) { plog.init(comm); serial = false; }
	void save(const std::string & fname)
	{
		if (serial)
			slog.save(fname);
		else
			plog.save(fname);
	}
protected:
	time_log<machine_mode::SERIAL> slog;
	time_log<machine_mode::MPI>    plog;
	bool active;
	bool serial;
};

}


#endif

