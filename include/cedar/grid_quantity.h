#ifndef CEDAR_GRID_QUANTITY_H
#define CEDAR_GRID_QUANTITY_H

#include <cedar/types.h>
#include <cedar/array_base.h>


namespace cedar {

/**
   Base class for objects that live on a grid.

   This provides useful iterators over the interior and ghost regions of a grid.
*/
template <typename len_type, unsigned short ND>
class grid_quantity : public virtual array_base<len_type>
{
public:

	/**
	   Provides the number of interior points (excluding ghosts) for a
	   given dimension.

	   @param i The dimension.
	   @returns Number of interior points for the requested dimension.
	*/
	virtual len_type shape(int i) const
	{
		return this->len(i) - 2*num_ghosts;
	}


	/**
	   Provides an iterator over the interior of the grid for a given
	   dimension.

	   @param i The dimension.
	   @returns Iterator over the grid's interior.
	*/
	const virtual range_t<len_type> & range(int i) const
	{
		#ifdef DEBUG
		return range_.at(i);
		#else
		return range_[i];
		#endif
	}


	/**
	   Provides an iterator over the entire grid (including ghosts)
	   for a given dimension.

	   @param i The dimension.
	   @returns Iterator over a grid dimension including ghost points.
	*/
	const virtual range_t<len_type> & grange(int i) const
	{
		#ifdef DEBUG
		return grange_.at(i);
		#else
		return grange_[i];
		#endif
	}


protected:
	int num_ghosts;
	std::array<range_t<len_type>, ND> range_;
	std::array<range_t<len_type>, ND> grange_;
};

}
#endif
