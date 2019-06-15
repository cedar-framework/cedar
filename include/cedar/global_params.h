#ifndef CEDAR_GLOBAL_PARAMS_H
#define CEDAR_GLOBAL_PARAMS_H

#include <memory>
#include <cstdint>
#include <iostream>

#include <cedar/config.h>

namespace cedar {

using loglevel_t = unsigned int;

namespace loglevel
{
constexpr loglevel_t status = 1;
constexpr loglevel_t info   = 1<<1;
constexpr loglevel_t error  = 1<<2;
constexpr loglevel_t memory = 1<<3;
constexpr loglevel_t debug  = 1<<4;
constexpr loglevel_t timer  = 1<<5;
}

enum class memtype { system, managed };

struct global_params
{
	loglevel_t log_level;
	memtype memory_type;
	friend std::ostream & operator<<(std::ostream & os, const global_params & obj);
};

loglevel_t getloglevel(config & conf);

global_params build_global_params(config & conf);

}

#endif
