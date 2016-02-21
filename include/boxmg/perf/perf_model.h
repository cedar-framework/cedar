#ifndef BOXMG_PERF_MODEL_H
#define BOXMG_PERF_MODEL_H

namespace boxmg {

class perf_model
{
public:
	virtual float time() = 0;
	virtual void set_comp_param(float tc) { this->tc = tc; }
	virtual void set_comm_param(float ts, float tw) { this->ts = ts; this->tw = tw; }

protected:
	float ts, tw, tc;
};

}

#endif
