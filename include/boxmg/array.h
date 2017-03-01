#ifndef BOXMG_ARRAY_H
#define BOXMG_ARRAY_H

#include <cassert>
#include <tuple>
#include <array>
#include <boxmg/types.h>
#include <boxmg/array_base.h>

namespace boxmg {

template <typename len_type, typename data_type, unsigned short ND>
class array : public virtual array_base<len_type>
{
protected:
	AlignedVector<data_type> vec;
	std::array<len_type, ND> strides;
	std::array<len_type, ND> extents;

public:

	len_type unpack_extents(len_type n)
	{
		extents[ND-1] = n;

		return ND-2;
	}

	template <typename... T> len_type unpack_extents(len_type n, T... args)
	{
		auto pos = unpack_extents(std::forward<decltype(args)>(args)...);
		extents[pos] = n;

		return pos-1;
	}

	array() {};
	template <typename... T> array(T... args)
	{
		reshape(std::forward<decltype(args)>(args)...);
	}


	template <typename... T> void reshape(T... args)
	{
		unpack_extents(std::forward<decltype(args)>(args)...);
		len_type len = 1;
		for (unsigned short i = 0; i < ND; i++)
			len *= extents[i];
		vec.resize(len);

		strides[0] = 1;
		for (unsigned short i = 1; i < ND; i++) {
			strides[i] = 1;
			for (unsigned short j = 0; j < i; j++) {
				strides[i] *= extents[j];
			}

		}
	}


	std::tuple<short, len_type> get_offset(len_type i) const
	{
		#ifdef BOUNDS_CHECK
		assert(i < extents[ND-1]);
		#endif
		return std::make_tuple(ND-2, i*strides[ND-1]);
	}


	template<typename... T> std::tuple<short, len_type> get_offset(len_type i, T... args) const
	{
		auto offset = get_offset(std::forward<decltype(args)>(args)...);
		auto pos = std::get<0>(offset);
		#ifdef BOUNDS_CHECK
		assert(pos >= 0);
		assert(i < extents[pos]);
		#endif
		return std::make_tuple(pos-1,
		                       std::get<1>(offset) + i*strides[pos]);
	}


	template<typename... T> data_type & operator()(T... args)
	{
		auto offset = get_offset(std::forward<decltype(args)>(args)...);
		#ifdef BOUNDS_CHECK
		assert(std::get<0>(offset) == -1);
		return vec.at(std::get<1>(offset));
		#else
		return vec[std::get<1>(offset)];
		#endif
	}


	template<typename... T> const data_type & operator()(T... args) const
	{
		auto offset = get_offset(std::forward<decltype(args)>(args)...);
		#ifdef BOUNDS_CHECK
		assert(std::get<0>(offset) == -1);
		return vec.at(std::get<1>(offset));
		#else
		return vec[std::get<1>(offset)];
		#endif
	}


	template<typename... T>	len_t index(T... args) const
	{
		auto offset = get_offset(std::forward<decltype(args)>(args)...);
		return std::get<1>(offset);
	}


	virtual len_type len(unsigned short i) const
	{
		#ifdef BOUNDS_CHECK
		assert(i < ND);
		#endif
		return extents[i];
	}


	virtual len_type & len(unsigned short i)
	{
		#ifdef BOUNDS_CHECK
		assert(i < ND);
		#endif
		return extents[i];
	}


	len_type stride(unsigned short i) const
	{
		#ifdef BOUNDS_CHECK
		assert(i < ND);
		#endif
		return strides[i];
	}


	len_type & stride(unsigned short i)
	{
		#ifdef BOUNDS_CHECK
		assert(i < ND);
		#endif
		return strides[i];
	}


	void set(data_type v)
	{
		for (auto&& val: vec)
			val = v;
	}


	void scale(data_type scalar)
	{
		for (auto&& val: vec)
			val *= scalar;
	}

	data_type * data() { return vec.data(); }
};

}
#endif
