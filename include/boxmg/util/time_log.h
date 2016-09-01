#ifndef BOXMG_TIMELOG_H
#define BOXMG_TIMELOG_H

#include <string>
#include <vector>

namespace boxmg {
	class time_log
	{
	public:
		time_log();
		void begin(std::string label);
		void end(std::string end);
		void up();
		void down();
		void add_comm(MPI_Comm comm) {comms.push_back(comm);}
		void redist(MPI_Comm comm) {add_comm(comm);redist_levels.push_back(lvl);}
		void save(std::string fname);

	private:
		std::vector<MPI_Comm> comms;
		std::vector<int> redist_levels;
		std::vector<std::map<std::string, double>> ltimes;
		std::vector<std::map<std::string, double>> stimes;
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
	inline void timer_redist(MPI_Comm comm) {tlog.redist(comm);}
	inline void timer_init(MPI_Comm comm){tlog.add_comm(comm);}
	inline void timer_save(std::string fname){tlog.save(fname);}
}

#endif
