#ifndef CEDAR_STENCIL_OP_ND_H
#define CEDAR_STENCIL_OP_ND_H

#include <iostream>
#include <utility>
#include <cedar/types.h>
#include <cedar/array.h>
#include <cedar/grid_quantity.h>

namespace cedar {

template <class T, T... Vs>
struct integer_sequence { };

template <class T, class, class, class = integer_sequence<T>, class = integer_sequence<T, 0>, class = void>
struct make_integer_sequence_impl;

template <class T, T ICV1, T... Res, T... Pow>
struct make_integer_sequence_impl<T, std::integral_constant<T, ICV1>, std::integral_constant<T, 0>, integer_sequence<T, Res...>, integer_sequence<T, Pow...>, typename std::enable_if<(ICV1 > 0)>::type>: make_integer_sequence_impl<T, std::integral_constant<T, ICV1/2>, std::integral_constant<T, ICV1%2>, integer_sequence<T, Res...>, integer_sequence<T, Pow..., (Pow + sizeof...(Pow))...>> { };

template <class T, T ICV1, T... Res, T... Pow>
struct make_integer_sequence_impl<T, std::integral_constant<T, ICV1>, std::integral_constant<T, 1>, integer_sequence<T, Res...>, integer_sequence<T, Pow...>, void>: make_integer_sequence_impl<T, std::integral_constant<T, ICV1/2>, std::integral_constant<T, ICV1%2>, integer_sequence<T, Pow..., (Res + sizeof...(Pow))...>, integer_sequence<T, Pow..., (Pow + sizeof...(Pow))...>> { };

template <class T, class Res, class Pow>
struct make_integer_sequence_impl<T, std::integral_constant<T, 0>, std::integral_constant<T, 0>, Res, Pow, void> {
   using type = Res;
};

template <class T, T V>
using make_integer_sequence = typename make_integer_sequence_impl<T, std::integral_constant<T, V/2>, std::integral_constant<T, V%2>>::type;

template <std::size_t V>
using make_index_sequence = make_integer_sequence<std::size_t, V>;

template <std::size_t... V>
using index_sequence = integer_sequence<std::size_t, V...>;

template<std::size_t nd, class stype, class = make_index_sequence<nd>>
	class stencil_op_nd;

template <std::size_t nd, typename stype, std::size_t... Is>
	class stencil_op_nd<nd, stype, index_sequence<Is...>> :
	public array<real_t, nd+1>,
	public grid_quantity<len_t, nd>
{
public:
stencil_op_nd() {}
stencil_op_nd(decltype(Is, len_t{})... args)
	{
		this->num_ghosts = 1;
		array<real_t, nd+1>::reshape(std::forward<decltype(args)>(args+2*this->num_ghosts)...,
		                                    static_cast<len_t>(stype::ndirs));

		for (std::size_t i = 0; i < nd; ++i) {
			this->range_[i] = cedar::range(static_cast<len_t>(this->num_ghosts),
			                               static_cast<len_t>(this->len(i)-1));
			this->grange_[i] = cedar::range(static_cast<len_t>(0),
			                                this->len(i));
		}
	}

	using array<real_t, nd+1>::len;
	using array<real_t, nd+1>::data;
	using array<real_t, nd+1>::set;
	using array<real_t, nd+1>::index;
	real_t & operator()(decltype(Is, len_t{})... args, stype dir)
	{
		return array<real_t,nd+1>::operator()(std::forward<decltype(args)>(args)...,static_cast<len_t>(dir));
	}


	real_t operator()(decltype(Is, len_t{})... args, stype dir) const
	{
		return array<real_t,nd+1>::operator()(std::forward<decltype(args)>(args)...,static_cast<len_t>(dir));
	}

	int ndirs() { return static_cast<int>(stype::ndirs); }
};

template <std::size_t nd, typename stype, std::size_t... Is>
	std::ostream & operator<<(std::ostream & os, const stencil_op_nd<nd,stype,index_sequence<Is...>> & obj);

template<class stype>
	struct stencil_ndirs
	{
		static const int value = static_cast<int>(stype::ndirs);
	};
}

#endif
