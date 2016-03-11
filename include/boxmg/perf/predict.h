#ifndef BOXMG_PERF_PREDICT_H
#define BOXMG_PERF_PREDICT_H

#include <vector>

namespace boxmg {
std::vector<int> predict_redist(int nprocx, int nprocy,
                                int ngx, int ngy);
}

#endif
