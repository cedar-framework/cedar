#ifndef BOXMG_CONST_MODEL_H
#define BOXMG_CONST_MODEL_H

#include <boxmg/types.h>
#include <boxmg/perf/perf_model.h>

namespace boxmg {

class const_model : public perf_model
{
public:
const_model(float t): tm(t) {}
	virtual float time() { return tm; }
	virtual void rep(std::ostream & os) const { os << "const model"; }

private:
	float tm;
};

}

#endif
