#ifndef BOXMG_PERF_MODEL_H
#define BOXMG_PERF_MODEL_H

namespace boxmg {

class perf_model
{
public:
	virtual float time() = 0;
protected:
	float ts, tw, tc;
};

}

#endif
