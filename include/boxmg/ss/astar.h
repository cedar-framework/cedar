#ifndef BOXMG_SS_ASTAR_H
#define BOXMG_SS_ASTAR_H

#include <queue>

#include <boxmg/ss/node.h>

namespace boxmg { namespace ss {

template<class solution_type, class problem_type>
solution_type astar(problem_type & problem)
{
	typedef typename problem_type::state_type state_type;
	typedef typename problem_type::action_type action_type;
	typedef typename problem_type::cost_type cost_type;
	using node_type = node<state_type,action_type,cost_type>;
	using node_ptr = std::shared_ptr<node_type>;

	auto init_node = std::make_shared<node_type>();
	init_node->state = problem.initial_state;
	init_node->path_cost = 0.0;

	auto heuristic = [](node_ptr nd) { return 0.0; };
	auto astar_priority = [heuristic](node_ptr lhs, node_ptr rhs) {
		cost_type lhs_cost = lhs->path_cost + heuristic(lhs);
		cost_type rhs_cost = rhs->path_cost + heuristic(rhs);

		return !(lhs_cost < rhs_cost);
	};
	std::priority_queue<node_ptr, std::vector<node_ptr>, decltype(astar_priority)> frontier(astar_priority);
	frontier.push(init_node);

	while (!frontier.empty()) {
		auto node = frontier.top(); frontier.pop();
		if (problem.goal_test(node->state)) return solution_type(node);
		// explored.insert(node.state);
		for (auto action : problem.actions(node->state)) {
			auto child = child_node(problem, node, action);
			// if (!(child.state in explored or child.state in frontier_map.keys())) {
			frontier.push(child);
			// else if child.state in frontier with higher cost {{
			// replace frontier node with child
		}
	}

	return solution_type(frontier.top());
}

}}

#endif
