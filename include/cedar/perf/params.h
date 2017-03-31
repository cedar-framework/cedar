#ifndef CEDAR_PERF_PARAMS_H
#define CEDAR_PERF_PARAMS_H

#include <cedar/config/reader.h>

namespace cedar {

class params
{
public:
	static float compute_tc(int nd, config::reader & conf);
};

}

#endif
