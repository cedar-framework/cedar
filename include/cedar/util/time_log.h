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
			if (stimes[lvl].find(label) == stimes[lvl].end()) {
				log::error << label << " timer never began!" << std::endl;
			}

			auto endtime = timer_util<mmode>::now();
			auto elapsed = timer_util<mmode>::duration(stimes[lvl][label], endtime);
			if (ltimes[lvl].find(label) == ltimes[lvl].end()) {
				ltimes[lvl][label] = elapsed;
				counts[lvl][label] = 1;
			} else {
				ltimes[lvl][label] += elapsed;
				counts[lvl][label] += 1;
			}
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

	template<> inline void time_log<machine_mode::MPI>::save(std::string fname)
	{
		using namespace boost::property_tree;
		ptree pt;
		const int TIMING = 0;
		const int ACTIVE_SIZE = 1;


		ptree children;
		int rank;
		MPI_Comm_rank(comm, &rank);

		{
			int max_levels_local = ltimes.size();
			int max_levels;
			MPI_Allreduce(&max_levels_local, &max_levels, 1, MPI_INT, MPI_MAX, comm);
			if (max_levels_local != max_levels)
				ltimes.resize(max_levels);
		}

		for (int i = ltimes.size() - 1; i >= 0; i--) {
			ptree child;
			int max_timers = ltimes[i].size();
			MPI_Bcast(&max_timers, 1, MPI_INT, 0, comm);
			std::vector<unsigned long int> hashes;
			hashes.resize(max_timers);
			int j = 0;
			for (auto &key : ltimes[i]) {
				hashes[j] = std::hash<std::string>{}(key.first);
				j++;
			}
			MPI_Bcast(hashes.data(), hashes.size(), MPI_UNSIGNED_LONG, 0, comm);
			for (auto &hash : hashes) {
				double max, min;
				double loc_time;
				double avgpkg[2];
				double loc_avgpkg[2];

				// Find local timing that matches this hash
				bool found = false;
				for (auto & timing : ltimes[i]) {
					if (hash == std::hash<std::string>{}(timing.first)) {
						double ratio, avg;
						loc_time = timing.second;
						MPI_Reduce(&loc_time, &min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
						MPI_Reduce(&loc_time, &max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
						loc_avgpkg[TIMING] = loc_time;
						loc_avgpkg[ACTIVE_SIZE] = 1;
						MPI_Reduce(loc_avgpkg, avgpkg, 2, MPI_DOUBLE, MPI_SUM, 0, comm);

						ratio = max / min;
						avg = avgpkg[TIMING] / avgpkg[ACTIVE_SIZE];

						if (rank == 0) {
							child.put(timing.first + ".min", min);
							child.put(timing.first + ".max", max);
							child.put(timing.first + ".ratio", ratio);
							child.put(timing.first + ".avg", avg);
							child.put(timing.first + ".count", counts[i][timing.first]);
						}

						found = true;
						break;
					}
				}
				if (not found) {
					loc_time = std::numeric_limits<double>::max();
					MPI_Reduce(&loc_time, &min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
					loc_time = std::numeric_limits<double>::min();
					MPI_Reduce(&loc_time, &max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
					loc_avgpkg[TIMING] = 0;
					loc_avgpkg[ACTIVE_SIZE] = 0;
					MPI_Reduce(loc_avgpkg, avgpkg, 2, MPI_DOUBLE, MPI_SUM, 0, comm);
				}
			}

			if (rank == 0)
				children.push_back(std::make_pair("", child));
		}

		if (rank == 0) {
			pt.add_child("levels", children);
			int size;
			MPI_Comm_size(comm, &size);
			pt.put("np", size);
			json_parser::write_json(fname, pt);
		}

		// Clean up (reinitialize timers)
		this->ltimes.clear();
		this->stimes.clear();
		this->ltimes.resize(1);
		this->stimes.resize(1);
		this->lvl = 0;
	}


	template<> inline void time_log<machine_mode::SERIAL>::save(std::string fname)
	{
		using namespace boost::property_tree;
		ptree pt;

		ptree children;

		for (int i = ltimes.size()-1; i >= 0; i--) {
			ptree child;
			for (auto &timing : ltimes[i]) {
				child.put(timing.first + ".time", timing.second);
				child.put(timing.first + ".count", counts[i][timing.first]);
			}
			children.push_back(std::make_pair("", child));
		}

		pt.add_child("levels", children);

		json_parser::write_json(fname, pt);

		// Clean up (reinitialize timers)
		this->ltimes.clear();
		this->stimes.clear();
		this->ltimes.resize(1);
		this->stimes.resize(1);
		this->lvl = 0;
	}

	extern time_log<machine_mode::MPI> tlog;
	extern time_log<machine_mode::SERIAL> tlog_ser;
	extern bool serial_timers;

	inline void timer_begin(std::string label) {
		if (serial_timers)
			tlog_ser.begin(label);
		else
			tlog.begin(label);
	}
	inline void timer_end(std::string label) {
		if (serial_timers)
			tlog_ser.end(label);
		else
			tlog.end(label);
	}
	inline void timer_up() {
		if (serial_timers)
			tlog_ser.up();
		else
			tlog.up();
	}
	inline void timer_down(){
		if (serial_timers)
			tlog_ser.down();
		else
			tlog.down();
	}
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
