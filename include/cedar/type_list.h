#ifndef CEDAR_TYPE_LIST_H
#define CEDAR_TYPE_LIST_H

namespace cedar
{
	template<class... T> struct type_list {};

	template<class ta, class tb>
	struct type_cat;

	template<typename ... a, typename ... b>
	struct type_cat<type_list<a ...>, type_list<b ...>>
	{ typedef type_list<a ..., b ...> type; };
}

#endif
