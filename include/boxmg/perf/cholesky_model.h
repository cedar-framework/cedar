#ifndef BOXMG_CHOLESKY_MODEL_H
#define BOXMG_CHOLESKY_MODEL_H

#include <boxmg/types.h>
#include <boxmg/perf/perf_model.h>

namespace boxmg {

class cholesky_model : public perf_model
{
public:
	cholesky_model(len_t n);
	virtual float time();
	virtual void rep(std::ostream & os) const;

private:
	len_t n;
};

}

#endif
