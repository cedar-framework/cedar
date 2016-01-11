#ifndef BOXMG_GRID_QUANTITY_H
#define BOXMG_GRID_QUANTITY_H

#include <boxmg/types.h>
#include <boxmg/array_base.h>


namespace boxmg {

template <typename len_type, unsigned short ND>
class grid_quantity : public virtual array_base<len_type>
{
public:
	virtual len_type shape(int i) const
	{
		return this->len(i) - 2*num_ghosts;
	}


	const virtual range_t<len_type> & range(int i) const
	{
		#ifdef DEBUG
		return range_.at(i);
		#else
		return range_[i];
		#endif
	}


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
