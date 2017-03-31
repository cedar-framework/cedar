#ifndef CEDAR_SS_SOLUTION_H
#define CEDAR_SS_SOLUTION_H

#include <cedar/ss/node.h>

namespace cedar { namespace ss {

template<class state_type, class action_type, class cost_type>
class solution
{
	typedef node<state_type,action_type,cost_type> node_type;
	using node_ptr = std::shared_ptr<node_type>;

public:
solution(node_ptr nd): sol_node(nd) {}

protected:
	node_ptr snode() { return sol_node; }
private:
	node_ptr sol_node;
};

}}

#endif
