#ifndef BOXMG_CORE_H
#define BOXMG_CORE_H

#include <boost/multi_array.hpp>
#include <vector>

namespace boxmg { namespace core {

template <typename D, unsigned short ND>
class array
{
private:
	boost::multi_array<D,ND> barr;

public:

	std::vector<int> unpack_extents()
	{
		std::vector<int> v;
		return v;
	}

	template <typename... T> std::vector<int> unpack_extents(int n, T... args)
	{
		auto v = unpack_extents(std::forward<decltype(args)>(args)...);
		v.push_back(n);
		return v;
	}

	array() {};
	template <typename... T> array(T... args)
		: barr(unpack_extents(std::forward<decltype(args)>(args)...))	{}

	template<typename IT> auto operator()(IT &&i) -> decltype((barr[i]))
	{
		return barr[i];
	}

	template<typename IT, typename... T>
		auto operator()(IT &&i, T... args) -> decltype(((*this)(args...)[i]))
	{
		return (*this)(std::forward<decltype(args)>(args)...)[i];
	}

	template <typename T> auto operator[](T&& i) -> decltype((barr[i]))
	{
		return barr[i];
	}


};

}}
#endif
