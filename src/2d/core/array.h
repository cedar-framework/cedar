#ifndef BOXMG_2D_CORE_ARRAY_H
#define BOXMG_2D_CORE_ARRAY_H

#include <iostream>
#include <cassert>
#include <stdlib.h>
#include <memory>
#include <array>

#include "boxmg-common.h"

namespace boxmg { namespace bmg2d {

template<typename S,typename D> class Array;
template<typename S, typename D> std::ostream & operator<< (std::ostream &os, const Array<S,D> & obj);

/**
 * Basic 2D array
 */
template <typename S, typename D>
class Array
{
public:
	Array() {};
	Array(S nx, S ny, unsigned int nghosts=1, bool create=true);

	D * data() { return data_.data(); };
	virtual D & operator()(S i, S j)
	{
		#ifdef DEBUG
		check_bounds(i,j);
		return data_.at(i*stride_[0]+j*stride_[1]);
		#else
		return data_.data()[i*stride_[0]+j*stride_[1]];
		#endif
	};

	const virtual D operator()(S i, S j) const
	{
		return data_.data()[i*stride_[0]+j*stride_[1]];
	};

	virtual S len(int i) const
	{
		#ifdef DEBUG
		return len_.at(i);
		#else
		return len_[i];
		#endif
	};

	virtual S shape(int i) const
	{
		#ifdef DEBUG
		return len_.at(i) - 2*nghosts_;
		#else
		return len_[i] - 2*nghosts_;
		#endif
	}

	virtual S stride(int i) const
	{
		#ifdef DEBUG
		return stride_.at(i);
		#else
		return stride_[i];
		#endif
	}

	virtual S & stride(int i)
	{
		#ifdef DEBUG
		return stride_.at(i);
		#else
		return stride_[i];
		#endif
	}

	inline const virtual S index(const S i, const S j) const { return i*stride_[0] + j*stride_[1]; };


	const virtual range_t<S> & range(int i) const
	{
		#ifdef DEBUG
		return range_.at(i);
		#else
		return range_[i];
		#endif
	}


	const virtual range_t<S> & grange(int i) const
	{
		#ifdef DEBUG
		return grange_.at(i);
		#else
		return grange_[i];
		#endif
	}


	virtual void check_bounds(S i, S j) const;

	static Array<S,D> ones(S nx, S ny)
	{
		Array<S,D> ret(nx,ny);

		for (auto j: ret.grange(1)) {
			for (auto i: ret.grange(0)) {
				ret(i,j) = 1.0;
			}
		}

		return ret;
	}


	static Array<S,D> zeros(S nx, S ny)
	{
		Array<S,D> ret(nx,ny);

		for (auto j: ret.grange(1)) {
			for (auto i: ret.grange(0)) {
				ret(i,j) = 0.0;
			}
		}

		return ret;
	}

	void scale(D v)
	{
		for (auto j: range_[1]) {
			for (auto i: range_[0]) {
				(*this)(i,j) = (*this)(i,j)*v;
			}
		}
	}

	void set(D v)
	{
		for (auto j: grange_[1]) {
			for (auto i: grange_[0]) {
				(*this)(i,j) = v;
			}
		}
	}

	AlignedVector<D> & vec() { return data_; }

	friend std::ostream& operator<< <> (std::ostream& os, const Array<S,D>& obj);

protected:
	AlignedVector<D> data_;
	std::array<S, 2> len_;
	std::array<S, 2> stride_;
	std::array<range_t<S>, 2> range_;
	std::array<range_t<S>, 2> grange_;
	unsigned int nghosts_;
};


template <typename S, typename D>
std::ostream & operator<< (std::ostream &os, const Array<S,D> & obj)
{
	for (S j = 0; j < obj.len(1); j++) {
		if (j == obj.nghosts_) os << '\n';
		for (S i = 0; i < obj.len(0); i++) {
			if (i == obj.nghosts_) os << ' ';
			os << std::to_string(obj(i,j)) << ' ';
			if (i == obj.len_[1] - 2*obj.nghosts_) os << ' ';
		}
		if (j == obj.len_[0] - 2*obj.nghosts_) os << '\n';
		os << '\n';
	}
	return os;
}

template <typename S, typename D>
	Array<S,D>::Array(S nx, S ny, unsigned int nghosts, bool create): nghosts_(nghosts)
{
	nx += 2*nghosts;
	ny += 2*nghosts;
	len_[0] = nx;
	len_[1] = ny;

	range_[0] = boxmg::range(static_cast<S>(nghosts), static_cast<S>(nx - nghosts));
	range_[1] = boxmg::range(static_cast<S>(nghosts), static_cast<S>(ny - nghosts));

	grange_[0] = boxmg::range(static_cast<S>(0), nx);
	grange_[1] = boxmg::range(static_cast<S>(0), ny);

	stride_[0] = 1;
	stride_[1] = nx;

	if (create) {
		data_.resize(len_[0]*len_[1]);
	}
}


template <typename S, typename D>
void Array<S,D>::check_bounds(S i, S j) const
{
	assert(i >=0 and j >=0 and i < len_[0] and j < len_[1]);
}

}}

#endif
