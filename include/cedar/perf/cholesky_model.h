#ifndef CEDAR_CHOLESKY_MODEL_H
#define CEDAR_CHOLESKY_MODEL_H

#include <cedar/types.h>
#include <cedar/perf/perf_model.h>

namespace cedar {

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
