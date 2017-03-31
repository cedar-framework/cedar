#ifndef CEDAR_TIMELOG_H
#define CEDAR_TIMELOG_H

#include <string>
#include <vector>
#include <map>
#include <cedar/mpi/redist_comms.h>

namespace cedar {
	class time_log
	{
	public:
		time_log();
		void begin(std::string label);
		void end(std::string end);
		void up();
		void down();
		void init(MPI_Comm comm) { this->comm = comm; }
		void redist(redist_comms comm) {rcomms.push_back(comm);redist_levels.push_back(lvl);}
		void save(std::string fname);

	private:
		MPI_Comm comm;
		std::vector<int> redist_levels;
		std::vector<redist_comms> rcomms;
		std::vector<std::map<std::string, double>> ltimes;
		std::vector<std::map<std::string, double>> stimes;
		std::vector<std::map<std::string, int>> counts;
		int lvl;
	};

	extern time_log tlog;

	inline void timer_begin(std::string label) {
		tlog.begin(label);
	}
	inline void timer_end(std::string label) {
		tlog.end(label);
	}
	inline void timer_up() {tlog.up();}
	inline void timer_down(){tlog.down();}
	inline void timer_redist(redist_comms comm) {tlog.redist(comm);}
	inline void timer_init(MPI_Comm comm){tlog.init(comm);}
	inline void timer_save(std::string fname){tlog.save(fname);}
}

#endif
