#ifndef CEDAR_ARRAY_H
#define CEDAR_ARRAY_H

#include <cassert>
#include <tuple>
#include <array>
#include <cedar/types.h>
#include <cedar/array_base.h>

#include <ftl/Buffer.hpp>
#include <cedar/device.h>

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
    ftl::FlatBuffer<data_type, len_type> base_buffer;

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

	aarray() {};

    // aarray(const aarray<len_type, data_type, ND>& other):
    //     base_buffer(other.base_buffer.get_numel()) {

    //     data_type* host_ptr = other.base_buffer.get_host_impl()->get_host_pointer();
    //     base_buffer.set_data(host_ptr, host_ptr + other.flat_len);

    //     reshape_array(other.extents);
    //     flat_len = other.flat_len;
    // }

    // aarray& operator=(const aarray<len_type, data_type, ND>& other) {
    //     base_buffer.resize(other.base_buffer.get_numel());
    //     data_type* host_ptr = other.base_buffer.get_host_impl()->get_host_pointer();
    //     base_buffer.set_data(host_ptr, host_ptr + other.flat_len);

    //     reshape_array(other.extents);
    //     flat_len = other.flat_len;

    //     return *this;
    // }

        template <typename... T> aarray(data_type *ext, T... args)
	{
		reshape(ext, std::forward<decltype(args)>(args)...);
	}

        template <typename... T> aarray(const ftl::BufferAllocateDevice device, data_type *ext, T... args):
        base_buffer(device)
	{
		reshape(ext, std::forward<decltype(args)>(args)...);
	}

	template <typename... T> aarray(T... args)
	{
		reshape(std::forward<decltype(args)>(args)...);
	}

        template <typename... T> aarray(const ftl::BufferAllocateDevice device, T... args):
        base_buffer(device)
	{
		reshape(std::forward<decltype(args)>(args)...);
	}

	template <typename... T> void init(T... args)
	{
		reshape(std::forward<decltype(args)>(args)...);
	}

    void reshape_array(const std::array<len_type, ND>& new_extents) {
		len_type len = 1;
		for (unsigned short i = 0; i < ND; i++) {
                    extents[i] = new_extents[i];
                    len *= extents[i];
                }
                base_buffer.resize(len);
		flat_len = len;
		base_ptr = base_buffer.data();

		strides[0] = 1;
		for (unsigned short i = 1; i < ND; i++) {
			strides[i] = 1;
			for (unsigned short j = 0; j < i; j++) {
				strides[i] *= extents[j];
			}

		}
	}

	template <typename... T> void reshape(T... args)
	{
		unpack_extents(std::forward<decltype(args)>(args)...);
		len_type len = 1;
		for (unsigned short i = 0; i < ND; i++)
			len *= extents[i];
                base_buffer.resize(len);
		flat_len = len;
		base_ptr = base_buffer.data();

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
		unpack_extents(std::forward<decltype(args)>(args)...);
		len_type len = 1;
		for (unsigned short i = 0; i < ND; i++)
			len *= extents[i];

		flat_len = len;
                base_buffer.set_ext_host_data(ext, flat_len);
                base_ptr = base_buffer.data();

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
		for (len_type i = 0; i < flat_len; i++)
			base_ptr[i] = v;

                base_buffer.mark_host_dirty();
	}


	void scale(data_type scalar)
	{
		for (len_type i = 0; i < flat_len; i++)
			base_ptr[i] *= scalar;

                base_buffer.mark_host_dirty();
	}

	data_type * data() {
            base_buffer.dev_to_host();
            return base_ptr;
        }

    ftl::FlatBuffer<data_type, len_type> to_flat_buffer() const {
        return base_buffer;
    }

    ftl::FlatBuffer<data_type, len_type>& to_flat_buffer() {
        return base_buffer;
    }

    ftl::Buffer<data_type, len_type> to_buffer() const {
        ftl::Buffer<data_type, len_type> buf(base_buffer);
        buf.reshape(std::vector<len_t>(extents.cbegin(), extents.cend()));
        return buf;
    }

    operator ftl::Buffer<data_type, len_type>() const {
        return to_buffer();
    }

    bool has_cpu() const {
        return base_buffer.has_cpu();
    }

    bool has_gpu() const {
        return base_buffer.has_gpu();
    }

    /**
     * Ensure that data is up-to-date on CPU, allocating buffers
     * and transferring data from GPU if necessary.
     */
    void ensure_cpu() {
        base_buffer.dev_to_host();
    }

    /**
     * Ensure that data is up-to-date on GPU, allocating buffers
     * and transferring data from GPU if necessary.
     */
    void ensure_gpu() {
        base_buffer.host_to_dev();
    }

    template <typename device>
    typename std::enable_if<std::is_same<device, cedar::cpu>::value>::type
    ensure() {
        ensure_cpu();
    }

    template <typename device>
    typename std::enable_if<std::is_same<device, cedar::gpu>::value>::type
    ensure() {
        ensure_gpu();
    }

    void mark_cpu_dirty(bool dirty=true) {
        base_buffer.mark_host_dirty(dirty);
    }

    void mark_gpu_dirty(bool dirty=true) {
        base_buffer.mark_device_dirty(dirty);
    }

    template <typename device>
    typename std::enable_if<std::is_same<device, cedar::cpu>::value>::type
    mark_dirty(bool dirty=true) {
        mark_cpu_dirty(dirty);
    }

    template <typename device>
    typename std::enable_if<std::is_same<device, cedar::gpu>::value>::type
    mark_dirty(bool dirty=true) {
        mark_gpu_dirty(dirty);
    }

    bool is_cpu_dirty() const {
        return base_buffer.is_host_dirty();
    }

    bool is_gpu_dirty() const {
        return base_buffer.is_dev_dirty();
    }
};

template<typename data_type, unsigned short ND>
	using array = aarray<len_t, data_type, ND>;

}
#endif
