#ifndef BOXMG_ARRAY_H
#define BOXMG_ARRAY_H

#include <cassert>
#include <tuple>

namespace boxmg {

template <typename D, unsigned short ND>
class array
{
private:
	D *arr;
	std::array<int, ND> strides;
	std::array<int, ND> extents;

public:

	std::vector<int> unpack_extents()
	{
		std::vector<int> v;
		return v;
	}

	int unpack_extents(int n)
	{
		extents[ND-1] = n;

		return ND-2;
	}

	template <typename... T> int unpack_extents(int n, T... args)
	{
		auto pos = unpack_extents(std::forward<decltype(args)>(args)...);
		extents[pos] = n;

		return pos-1;
	}

	array() {};
	template <typename... T> array(T... args)
	{
		auto pos = unpack_extents(std::forward<decltype(args)>(args)...);
		#ifdef BOUNDS_CHECK
		assert(pos == -1);
		#endif
		int len = 1;
		for (int i = 0; i < ND; i++)
			len *= extents[i];
		arr = new D[len];

		strides[0] = 1;
		for (int i = 1; i < ND; i++) {
			strides[i] = 1;
			for (int j = 0; j < i; j++) {
				strides[i] *= extents[j];
			}

		}
	}


	~array()
	{
		delete[] arr;
	}


	std::tuple<int, D*> get_base(int i)
	{
		#ifdef BOUNDS_CHECK
		assert(i < extents[ND-1]);
		#endif
		return std::make_tuple(ND-2, arr + i*strides[ND-1]);
	}


	template<typename... T> std::tuple<int, D*> get_base(int i, T... args)
	{
		auto base = get_base(std::forward<decltype(args)>(args)...);
		auto pos = std::get<0>(base);
		#ifdef BOUNDS_CHECK
		assert(pos >= 0);
		assert(i < extents[pos]);
		#endif
		return std::make_tuple(pos-1,
		                       std::get<1>(base) + i*strides[pos]);
	}

	template<typename... T> D & operator()(T... args)
	{
		auto base = get_base(std::forward<decltype(args)>(args)...);
		#ifdef BOUNDS_CHECK
		assert(std::get<0>(base) == -1);
		#endif
		D * ret = std::get<1>(base);
		return *ret;
	}

	int len(int i) const
	{
		#ifdef BOUNDS_CHECK
		assert(i < ND);
		#endif
		return extents[i];
	}

	int stride(int i) const
	{
		#ifdef BOUNDS_CHECK
		assert(i < ND);
		#endif
		return strides[i];
	}
};

}
#endif
