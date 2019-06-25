#ifndef CEDAR_ARRAY_H
#define CEDAR_ARRAY_H

#include <cassert>
#include <tuple>
#include <array>
#include <cedar/memory.h>
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
	bool owns_data;
	data_type *base_ptr;
	len_type flat_len;
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

	aarray() : owns_data(false) {};
	~aarray() { if (owns_data) memory::free(base_ptr); }
	aarray(const aarray &other)
		: owns_data(true), flat_len(other.flat_len),
		  strides(other.strides), extents(other.extents)
	{
		base_ptr = memory::alloc<data_type>(other.flat_len);
		std::copy(other.base_ptr, other.base_ptr + other.flat_len, base_ptr);
	}
	aarray(aarray&& other) noexcept {
		*this = std::move(other);
	}
	aarray & operator=(const aarray & other)
	{
		return *this = aarray(other);
	}
	aarray & operator=(aarray&& other) noexcept {
		if (this == &other) return *this;
		strides = std::move(other.strides);
		extents = std::move(other.extents);
		flat_len = other.flat_len;
		base_ptr = other.base_ptr;
		owns_data = other.owns_data;
		other.owns_data = false;
		other.base_ptr = nullptr;
		other.flat_len = 0;
		return *this;
	}
	template <typename... T> aarray(data_type *ext, T... args)
	{
		reshape(ext, std::forward<decltype(args)>(args)...);
	}

	template <typename... T> aarray(T... args)
	{
		reshape(std::forward<decltype(args)>(args)...);
	}

	template <typename... T> void init(T... args)
	{
		reshape(std::forward<decltype(args)>(args)...);
	}

	template <typename... T> void reshape(T... args)
	{
		unpack_extents(std::forward<decltype(args)>(args)...);
		len_type len = 1;
		for (unsigned short i = 0; i < ND; i++)
			len *= extents[i];

		base_ptr = memory::alloc<data_type>(len);
		owns_data = true;
		flat_len = len;

		strides[0] = 1;
		for (unsigned short i = 1; i < ND; i++) {
			strides[i] = 1;
			for (unsigned short j = 0; j < i; j++) {
				strides[i] *= extents[j];
			}

		}
	}


	template <typename... T> void reshape(data_type *ext, T... args)
	{
		base_ptr = ext;
		unpack_extents(std::forward<decltype(args)>(args)...);
		len_type len = 1;
		for (unsigned short i = 0; i < ND; i++)
			len *= extents[i];

		flat_len = len;
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
		return base_ptr[offset];
		#else
		return base_ptr[offset];
		#endif
	}


	template<typename... T> const data_type & operator()(T... args) const
	{
		len_type offset; short pos;
		std::tie(pos, offset) = get_offset(std::forward<decltype(args)>(args)...);
		#ifdef BOUNDS_CHECK
		assert(pos == -1);
		return base_ptr[offset];
		#else
		return base_ptr[offset];
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
		#ifdef OFFLOAD
		bool ongpu = memory::hint(base_ptr) == memory::location::gpu;
		#pragma omp target teams distribute parallel for simd if (ongpu)
		for (len_type i = 0; i < flat_len; i++)
			base_ptr[i] = v;
		#else
		for (len_type i = 0; i < flat_len; i++)
			base_ptr[i] = v;
		#endif
	}


	void scale(data_type scalar)
	{
		for (len_type i = 0; i < flat_len; i++)
			base_ptr[i] *= scalar;
	}

	data_type * data() { return base_ptr; }
	const data_type * cdata() const { return base_ptr; }
	len_type size() const { return flat_len; }
};

template<typename data_type, unsigned short ND>
	using array = aarray<len_t, data_type, ND>;

}
#endif
