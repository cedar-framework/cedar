#ifndef CEDAR_SS_PROBLEM_H
#define CEDAR_SS_PROBLEM_H

#include <memory>
#include <vector>

namespace cedar { namespace ss {

template<class T, class U, class V>
struct problem
{
	typedef T state_type;
	typedef U action_type;
	typedef V cost_type;

	virtual bool goal_test(state_type & state) = 0;
	virtual state_type result(state_type & state, action_type action) = 0;
	virtual cost_type step_cost(state_type & state, action_type action) = 0;
	virtual std::vector<action_type> actions(state_type & state) = 0;
	state_type initial_state;
};

}}
#endif
