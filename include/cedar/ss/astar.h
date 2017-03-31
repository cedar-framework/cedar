#ifndef CEDAR_SS_ASTAR_H
#define CEDAR_SS_ASTAR_H

#include <queue>

#include <cedar/ss/node.h>

namespace cedar { namespace ss {

template<class solution_type, class problem_type>
solution_type astar(problem_type & problem,
                    std::function<typename problem_type::cost_type(std::shared_ptr<node<typename problem_type::state_type, typename problem_type::action_type, typename problem_type::cost_type>>)> heuristic)
{
	typedef typename problem_type::state_type state_type;
	typedef typename problem_type::action_type action_type;
	typedef typename problem_type::cost_type cost_type;
	using node_type = node<state_type,action_type,cost_type>;
	using node_ptr = std::shared_ptr<node_type>;

	auto init_node = std::make_shared<node_type>();
	init_node->state = problem.initial_state;
	init_node->path_cost = 0.0;

	auto astar_priority = [heuristic](node_ptr lhs, node_ptr rhs) {
		cost_type lhs_cost = lhs->path_cost + heuristic(lhs);
		cost_type rhs_cost = rhs->path_cost + heuristic(rhs);

		return !(lhs_cost < rhs_cost);
	};
	std::priority_queue<node_ptr, std::vector<node_ptr>, decltype(astar_priority)> frontier(astar_priority);
	std::map<state_type, cost_type> explored;
	frontier.push(init_node);
	explored[init_node->state] = 0;

	while (!frontier.empty()) {
		node_ptr node;
		do {
			node = frontier.top(); frontier.pop();
		} while (node->path_cost > explored[node->state]);
		if (problem.goal_test(node->state)) return solution_type(node);
		for (auto action : problem.actions(node->state)) {
			auto child = child_node(problem, node, action);
			auto child_cost = child->path_cost;
			if ((explored.find(child->state) == explored.end()) or (child_cost < explored[child->state])) {
				explored[child->state] = child_cost;
				frontier.push(child);
			}
		}
	}

	return solution_type(frontier.top());
}

}}

#endif
