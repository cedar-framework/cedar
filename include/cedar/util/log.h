#ifndef CEDAR_LOG_H
#define CEDAR_LOG_H

#include <mpi.h>
#include <memory>
#include <map>
#include <string>
#include <iostream>
#include <ctime>
#include <sstream>

#include <iomanip>

#include "color_mod.h"
#include <cedar/config.h>

namespace cedar { namespace log {

void init();

unsigned int & lvl();
void set_comm(MPI_Comm comm);
MPI_Comm get_comm();
void push_comm(MPI_Comm comm);
void pop_comm();
void set_header_msg(std::string header_msg);
void init_level(config & conf);
void init(config & conf);
void push_level(std::string header, config & conf);
void pop_level();

std::string header();

class LogLevelBuf : public std::stringbuf
{
public:
	template<class T, class D>
	LogLevelBuf(T&& name, D&& colormod): name(std::forward<T>(name)),
		color(std::forward<D>(colormod)),
		ost(std::cout){}
	template<class T, class D>
    LogLevelBuf(T&& name, D&& colormod, std::ostream & ost): name(std::forward<T>(name)),
		color(std::forward<D>(colormod)),
		ost(ost) {}
	~LogLevelBuf() { pubsync(); }
	int sync();

private:
	std::string name;
	Color::Modifier color;
	std::ostream & ost;
};

class LevelLogger : public std::ostream
{
public:
	template<class T, class C>
		LevelLogger(T&& name, C&& colormod) : std::ostream(new LogLevelBuf(name,
		                                                                   std::forward<C>(colormod))), name(name){}
	template<class T, class C>
	LevelLogger(T&& name, C&& colormod, std::ostream &ost):
	std::ostream(new LogLevelBuf(std::forward<T>(name),
	                             std::forward<C>(colormod),
	                             ost)) {}
	~LevelLogger() { delete rdbuf(); }
	bool active();

private:
	std::string name;
};

extern LevelLogger memory;
extern LevelLogger error;
extern LevelLogger info;
extern LevelLogger status;
extern LevelLogger debug;
extern LevelLogger timer;

}}
#endif
