#ifndef CEDAR_CXX14_H
#define CEDAR_CXX14_H

#include <memory>
#include <tuple>

namespace std {
template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
	return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

template<class T, std::size_t N, class... Args>
struct tuple_id_from_type
{
	static constexpr auto value = N;
};

template<class T, std::size_t N, class... Args>
struct tuple_id_from_type<T, N, T, Args...>
{
	static constexpr auto value = N;
};

template<class T, std::size_t N, class U, class... Args>
struct tuple_id_from_type<T, N, U, Args...>
{
	static constexpr auto value = tuple_id_from_type<T, N+1, Args...>::value;
};

template<class T, class... Types>
T& get(tuple<Types...>& t) noexcept
{
	return get<tuple_id_from_type<T, 0, Types...>::value>(t);
}

template<class T, class... Types>
T&& get(tuple<Types...>&& t) noexcept
{
	std::forward<T&&>(get<tuple_id_from_type<T, 0, Types...>::value>(t));
}
}

#endif
