#include <boxmg/perf/cholesky_model.h>

using namespace boxmg;

cholesky_model::cholesky_model(len_t n): n(n) {}


float cholesky_model::time() const
{
	return (n*n)*tc;
}


void cholesky_model::rep(std::ostream &os) const
{
	os << "======== Cholesky model ========" << '\n';
	os << "size: " << n << '\n';
	os << "time: " << time() << '\n';
}
