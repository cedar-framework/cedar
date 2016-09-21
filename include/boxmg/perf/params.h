#ifndef BOXMG_PERF_PARAMS_H
#define BOXMG_PERF_PARAMS_H

#include <boxmg/config/reader.h>

namespace boxmg {

class params
{
public:
	static float compute_tc(int nd, config::reader & conf);
};

}

#endif
