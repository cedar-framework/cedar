#ifndef CEDAR_UTIL_BASIC_LOGGER_H
#define CEDAR_UTIL_BASIC_LOGGER_H

#include <mpi.h>
#include <iostream>
#include <stack>
#include <tuple>

#include <cedar/global_params.h>
#include "color_mod.h"

namespace cedar {


class LogLevelBuf : public std::stringbuf
{
public:
	template<class D>
	LogLevelBuf(loglevel_t levelid, loglevel_t & curlevel, std::string &header_msg,
	            MPI_Comm &comm, D&& colormod):
		levelid(levelid), curlevel(curlevel),
		header_msg(header_msg), comm(comm),
		color(std::forward<D>(colormod)),
		ost(std::cout){}
	template<class D>
	LogLevelBuf(loglevel_t levelid, loglevel_t & curlevel,
	            std::string & header_msg, MPI_Comm &comm,
	            D&& colormod, std::ostream & ost):
		levelid(levelid), curlevel(curlevel),
		header_msg(header_msg), comm(comm),
		color(std::forward<D>(colormod)),
		ost(ost) {}
	~LogLevelBuf() { pubsync(); }
	int sync();
	std::string header();

private:
	loglevel_t levelid;
	loglevel_t & curlevel;
	std::string & header_msg;
	MPI_Comm & comm;
	Color::Modifier color;
	std::ostream & ost;
};


class levellogger : public std::ostream
{
public:
	template<class C>
	levellogger(loglevel_t levelid, loglevel_t & curlevel,
	            std::string &header_msg, MPI_Comm &comm, C&& colormod) :
		std::ostream(new LogLevelBuf(levelid, curlevel, header_msg, comm,
		                             std::forward<C>(colormod))), levelid(levelid), curlevel(curlevel) {}
	template<class C>
	levellogger(loglevel_t levelid, loglevel_t & curlevel,
	            std::string &header_msg, MPI_Comm & comm,
	            C&& colormod, std::ostream &ost):
		std::ostream(new LogLevelBuf(levelid, curlevel,
		                             header_msg, comm,
		                             std::forward<C>(colormod),
		                             ost)),
		levelid(levelid), curlevel(curlevel) {}
	~levellogger() { delete rdbuf(); }
	bool active() { return (levelid & curlevel); }

protected:
	loglevel_t levelid;
	loglevel_t & curlevel;
};


class basic_logger
{
protected:
	loglevel_t level;
	std::string header_msg;
	MPI_Comm comm;
	std::stack<std::tuple<loglevel_t, std::string>> saved_levels;
	std::stack<MPI_Comm> saved_comms;
public:
	basic_logger(global_params & params);
	void set_header(std::string msg) { header_msg = msg; }
	void set_comm(MPI_Comm new_comm) { comm = new_comm; }
	loglevel_t & lvl() { return level; }
	MPI_Comm get_comm() { return comm; }
	void push_comm(MPI_Comm comm);
	void pop_comm();
	void push_level(std::string header, loglevel_t level);
	void pop_level();

	levellogger status;
	levellogger error;
	levellogger info;
	levellogger debug;
	levellogger timer;
	levellogger memory;
};

}

#endif
