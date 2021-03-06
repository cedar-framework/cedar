#ifndef CEDAR_PERF_MODEL_H
#define CEDAR_PERF_MODEL_H

#include <iostream>
#include <boost/property_tree/ptree.hpp>


namespace cedar {

class perf_model
{
public:
	virtual float time() const = 0;
	virtual void set_comp_param(float tc) { this->tc = tc; }
	virtual float get_comp_param() { return this->tc; }
	virtual void set_comm_param(float ts, float tw) { this->ts = ts; this->tw = tw; }
	virtual void rep(std::ostream &os) const = 0;
	virtual int nproc() const = 0;
	virtual void recur_times(boost::property_tree::ptree &) const = 0;
	friend std::ostream & operator<<(std::ostream &os, const perf_model & pm);


protected:
	float ts, tw, tc;
};

inline std::ostream & operator<<(std::ostream &os, const perf_model & pm) {
	pm.rep(os);
	return os;
}

}

#endif
