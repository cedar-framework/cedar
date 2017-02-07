#ifndef BOXMG_CHOLESKY_MODEL_H
#define BOXMG_CHOLESKY_MODEL_H

#include <boxmg/types.h>
#include <boxmg/perf/perf_model.h>

namespace boxmg {

class cholesky_model : public perf_model
{
public:
	cholesky_model(len_t n);
	virtual float time() const;
	virtual void rep(std::ostream & os) const;
	virtual void recur_times(boost::property_tree::ptree &) const {};
	virtual int nproc() const { return 1; }

private:
	len_t n;
};

}

#endif
