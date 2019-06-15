#ifndef CEDAR_LOG_H
#define CEDAR_LOG_H

#include <ostream>
#include <cedar/global_manager.h>

namespace cedar { namespace log {

using logger = global_manager<reg_globals>::logger;

inline void set_comm(MPI_Comm comm) { gman.get<logger>().set_comm(comm); }
inline MPI_Comm get_comm() { return gman.get<logger>().get_comm(); }
inline void push_comm(MPI_Comm comm) { gman.get<logger>().push_comm(comm); }
inline void pop_comm() { gman.get<logger>().pop_comm(); }
inline void set_header_msg(std::string header_msg) { gman.get<logger>().set_header(header_msg); }
inline void init_level(config & conf) {	gman.get<logger>().lvl() = getloglevel(conf); }
inline void push_level(std::string header, config & conf)
{
	gman.get<logger>().push_level(header, getloglevel(conf));
}
inline void pop_level() { gman.get<logger>().pop_level(); }
inline loglevel_t & lvl() { return gman.get<logger>().lvl(); }

class logwrapper
{
public:
	using wraptype = decltype(gman.get<logger>().status);

	logwrapper() : logwrap(nullptr) {}

	void wrap(wraptype *val)
	{
		logwrap = val;
	}

	bool active()
	{
		if (not logwrap) return false;
		else return logwrap->active();
	}

	template<typename T>
	logwrapper & operator<<(T&& x)
	{
		if (logwrap)
			*logwrap << std::forward<T>(x);
		return *this;
	}
	logwrapper & operator<<(std::ostream& (*manip)(std::ostream&))
	{
		if (logwrap)
			*logwrap << manip;
		return *this;
	}

protected:
	wraptype *logwrap;
};

extern logwrapper status;
extern logwrapper error;
extern logwrapper info;
extern logwrapper debug;
extern logwrapper timer;
extern logwrapper memory;

}}

#endif

