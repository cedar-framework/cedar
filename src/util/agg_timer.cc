#include <cedar/util/agg_timer.h>

namespace cedar
{

agg_timer::agg_timer(global_params & params) : active(true), serial(true) {}


template<> void time_log<machine_mode::MPI>::save(const std::string & fname)
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


template<> void time_log<machine_mode::SERIAL>::save(const std::string & fname)
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


}
