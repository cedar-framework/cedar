#ifndef CEDAR_ARRAY_H
#define CEDAR_ARRAY_H

#include <cassert>
#include <tuple>
#include <array>
#include <cedar/types.h>
#include <cedar/array_base.h>

namespace cedar {

/**
   Multidimensional array.

   Provides fortran-like ND array.  This is the core data structure
   used to store grid functions and stencil operators.  The inner
   index is the fastest moving (Fortran ordering).

   @tparam len_type Type for the lengths of the array's dimensions.
   @tparam data_type Type for what is stored in the array.
   @tparam ND Number of dimensions.
 */
template <typename len_type, typename data_type, unsigned short ND>
class aarray : public virtual array_base<len_type>
{
protected:
	AlignedVector<data_type> vec;     /** Vector where the data is stored. */
	std::array<len_type, ND> strides; /** Strides of each dimension, e.g., first index is stride 1. */
	std::array<len_type, ND> extents; /** Lengths of each dimension. */

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

	aarray() {};
	template <typename... T> aarray(T... args)
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


	std::pair<short, len_type> get_offset(len_type i) const
	{
		#ifdef BOUNDS_CHECK
		assert(i < extents[ND-1]);
		#endif
		return std::make_pair(ND-2, i*strides[ND-1]);
	}


	template<typename... T> std::pair<short, len_type> get_offset(len_type i, T... args) const
	{
		len_type offset; short pos;
		std::tie(pos, offset) = get_offset(std::forward<decltype(args)>(args)...);
		#ifdef BOUNDS_CHECK
		assert(pos >= 0);
		assert(i < extents[pos]);
		#endif
		return std::make_pair(pos-1,
		                      offset + i*strides[pos]);
	}


	template<typename... T> data_type & operator()(T... args)
	{
		len_type offset; short pos;
		std::tie(pos, offset) = get_offset(std::forward<decltype(args)>(args)...);
		#ifdef BOUNDS_CHECK
		assert(pos == -1);
		return vec.at(offset);
		#else
		return vec[offset];
		#endif
	}


	template<typename... T> const data_type & operator()(T... args) const
	{
		len_type offset; short pos;
		std::tie(pos, offset) = get_offset(std::forward<decltype(args)>(args)...);
		#ifdef BOUNDS_CHECK
		assert(pos == -1);
		return vec.at(offset);
		#else
		return vec[offset];
		#endif
	}


	template<typename... T>	len_t index(T... args) const
	{
		len_t offset;
		std::tie(std::ignore, offset) = get_offset(std::forward<decltype(args)>(args)...);
		return offset;
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

template<typename data_type, unsigned short ND>
	using array = aarray<len_t, data_type, ND>;

}
#endif
