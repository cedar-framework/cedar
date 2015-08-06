#ifndef BOXMG_UTIL_TYPES_H
#define BOXMG_UTIL_TYPES_H

#include <memory>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/foreach.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <cassert>

#include "align_allocator.h"

namespace std {
	template<typename T, typename ...Args>
		std::unique_ptr<T> make_unique( Args&& ...args )
	{
		return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
	}
}

namespace boxmg {

using len_t = unsigned int;
using real_t = double;

template <class T>
	using AlignedVector = std::vector<T, AlignAllocator<T,16>>;


template <class T>
boost::iterator_range<boost::counting_iterator<T> > range(T to)
{
    //these assertions are somewhat problematic:
    //might produce warnings, if T is unsigned
    assert(T() <= to);
    return boost::make_iterator_range(boost::counting_iterator<T>(0), boost::counting_iterator<T>(to));
}

template <class T>
boost::iterator_range<boost::counting_iterator<T> > range(T from, T to)
{
    assert(from <= to);
    return boost::make_iterator_range(boost::counting_iterator<T>(from), boost::counting_iterator<T>(to));
}

//iterator that can do increments in steps (positive and negative)
template <class T>
class range_iterator:
    public boost::iterator_facade<range_iterator<T>, const T, std::forward_iterator_tag>
{
    T value, incr;
public:
    range_iterator(T value, T incr = T()): value(value), incr(incr) {}
private:
    friend class boost::iterator_core_access;
    void increment() { value += incr; }
    bool equal(const range_iterator& other) const
    {
        //this is probably somewhat problematic, assuming that the "end iterator"
        //is always the right-hand value?
        return (incr >= 0 && value >= other.value) || (incr < 0 && value <= other.value);
    }
    const T& dereference() const { return value; }
};

template <class T>
boost::iterator_range<range_iterator<T> > range(T from, T to, T increment)
{
    assert((increment >= T() && from <= to) || (increment < T() && from >= to));
    return boost::make_iterator_range(range_iterator<T>(from, increment), range_iterator<T>(to));
}

template <class T>
using range_t = boost::iterator_range<boost::counting_iterator<T>>;
}

#endif
