#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <boxmg/util/log.h>
#include <boxmg/util/time_log.h>
#include <iostream>

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

	ptree children;
	for (int i = ltimes.size()-1; i >= 0; i--) {
		ptree child;
		for (auto &timing : ltimes[i]) {
			int size;
			double max, min, ratio, avg, tot;
			double loc_time;
			MPI_Comm comm = comms[0];
			loc_time = timing.second;

			MPI_Comm_size(comm, &size);

			MPI_Reduce(&loc_time, &min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
			MPI_Reduce(&loc_time, &max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
			MPI_Reduce(&loc_time, &tot, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

			ratio = max / min;
			avg = tot / size;

			child.put(timing.first + ".min", min);
			child.put(timing.first + ".max", max);
			child.put(timing.first + ".ratio", ratio);
			child.put(timing.first + ".avg", avg);

		}
		children.push_back(std::make_pair("", child));
	}
	pt.add_child("levels", children);

	int rank;
	MPI_Comm_rank(comms[0], &rank);
	if (rank == 0)
		json_parser::write_json(fname, pt);
}
