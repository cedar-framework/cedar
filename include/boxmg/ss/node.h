#ifndef BOXMG_SS_NODE_H
#define BOXMG_SS_NODE_H

#include <memory>

#include <iostream>
#include <boxmg/ss/problem.h>

namespace boxmg { namespace ss {

template<class T, class U, class V>
struct node
{
	typedef T state_type;
	typedef U action_type;
	typedef V cost_type;

node() : parent(NULL) {}

	state_type state;
	std::shared_ptr<node> parent;
	action_type action;
	cost_type path_cost;
};


template<class state_type, class action_type, class cost_type>
	std::shared_ptr<node<state_type,action_type,cost_type>> child_node(
		problem<state_type,action_type,cost_type> & problem,
	std::shared_ptr<node<state_type,action_type,cost_type>> pnode,
	action_type action)
{
	using node_type = node<state_type,action_type,cost_type>;
	auto ret = std::make_shared<node_type>();

	ret->state = problem.result(pnode->state, action);
	ret->parent = pnode;
	ret->action = action;
	ret->path_cost = pnode->path_cost + problem.step_cost(pnode->state, action);

	return ret;
}

}}


#endif
