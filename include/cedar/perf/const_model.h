#ifndef CEDAR_CONST_MODEL_H
#define CEDAR_CONST_MODEL_H

#include <cedar/types.h>
#include <cedar/perf/perf_model.h>

namespace cedar {

class const_model : public perf_model
{
public:
const_model(float t): tm(t) {}
	virtual float time() const { return tm; }
	virtual int nproc() const { return 1; }
	virtual void rep(std::ostream & os) const { os << "const model"; }
	virtual void recur_times(boost::property_tree::ptree &) const {};

private:
	float tm;
};

}

#endif
