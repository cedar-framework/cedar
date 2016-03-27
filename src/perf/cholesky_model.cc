#include <boxmg/perf/cholesky_model.h>

using namespace boxmg;

cholesky_model::cholesky_model(len_t n): n(n) {}


float cholesky_model::time()
{
	return (n*n)*tc;
}
