#include <iostream>
#include <limits>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <boxmg/util/log.h>
#include <boxmg/util/time_log.h>

using namespace boxmg;

namespace boxmg {
	time_log tlog;
}


time_log::time_log() : lvl(0)
{
	stimes.resize(1);
	ltimes.resize(1);
}


void time_log::up()
{
	lvl--;
}


void time_log::down()
{
	lvl++;
	if (lvl+1 > stimes.size()) {
		stimes.resize(stimes.size()+1);
		ltimes.resize(ltimes.size()+1);
	}
}


void time_log::begin(std::string label)
{
	stimes[lvl][label] = MPI_Wtime();
}


void time_log::end(std::string label)
{
	if (stimes[lvl].find(label) == stimes[lvl].end()) {
		log::error << label << " timer never began!" << std::endl;
	}

	double endtime = MPI_Wtime();
	double elapsed = endtime - stimes[lvl][label];
	if (ltimes[lvl].find(label) == ltimes[lvl].end()) {
		ltimes[lvl][label] = elapsed;
	} else {
		ltimes[lvl][label] += elapsed;
	}


	stimes[lvl][label] = 0;
}


void time_log::save(std::string fname)
{
	using namespace boost::property_tree;
	ptree pt;
	const int LEVEL_COUNT = 0;
	const int TIMING_COUNT = 1;
	const int TIMING = 0;
	const int ACTIVE_SIZE = 1;


	ptree children;

	int linfo[2];
	int loc_linfo[2];
	loc_linfo[LEVEL_COUNT] = ltimes.size();
	loc_linfo[TIMING_COUNT] = 0;
	for (auto i = 0; i < ltimes.size(); ++i) {
		loc_linfo[TIMING_COUNT] += ltimes[i].size();
	}
	MPI_Allreduce(&loc_linfo, &linfo, 2, MPI_INT, MPI_MAX, comm);

	for (auto i = 0; i < linfo[TIMING_COUNT] - loc_linfo[TIMING_COUNT]; ++i) {
		double max, min, ratio, avg;
		double loc_time;
		double avgpkg[2];
		double loc_avgpkg[2];

		loc_time = std::numeric_limits<double>::max();
		MPI_Reduce(&loc_time, &min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
		loc_time = std::numeric_limits<double>::min();
		MPI_Reduce(&loc_time, &max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
		loc_avgpkg[TIMING] = 0;
		loc_avgpkg[ACTIVE_SIZE] = 0;
		MPI_Reduce(loc_avgpkg, avgpkg, 2, MPI_DOUBLE, MPI_SUM, 0, comm);
	}

	for (int i = ltimes.size()-1; i >= 0; i--) {
		ptree child;
		for (auto &timing : ltimes[i]) {
			int size;
			double max, min, ratio, avg;
			double loc_time;
			double avgpkg[2];
			double loc_avgpkg[2];
			loc_time = timing.second;

			MPI_Comm_size(comm, &size);

			MPI_Reduce(&loc_time, &min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
			MPI_Reduce(&loc_time, &max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
			loc_avgpkg[TIMING] = loc_time;
			loc_avgpkg[ACTIVE_SIZE] = 1;
			MPI_Reduce(loc_avgpkg, avgpkg, 2, MPI_DOUBLE, MPI_SUM, 0, comm);

			ratio = max / min;
			avg = avgpkg[TIMING] / avgpkg[ACTIVE_SIZE];

			child.put(timing.first + ".min", min);
			child.put(timing.first + ".max", max);
			child.put(timing.first + ".ratio", ratio);
			child.put(timing.first + ".avg", avg);

		}
		children.push_back(std::make_pair("", child));
	}
	pt.add_child("levels", children);

	int rank, size;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);
	pt.put("np", size);
	if (rank == 0)
		json_parser::write_json(fname, pt);
}
