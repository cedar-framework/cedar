#ifndef BOXMG_2D_CORE_RELAX_STENCIL_H
#define BOXMG_2D_CORE_RELAX_STENCIL_H

#include "core/types.h"
#include "core/array.h"

namespace boxmg { namespace bmg2d { namespace core {


class RelaxStencil : public Array<len_t, real_t>
{
public:
	RelaxStencil() {};
	RelaxStencil(len_t nx, len_t ny, unsigned int nghosts=1);

	using Array<len_t,real_t>::index;

	len_t index(len_t i, len_t j, int dir) const
	{
		return dir*len_[0]*len_[1] + i*stride_[0] + j * stride_[1];
	}

	real_t operator()(len_t i, len_t j, int dir) const
	{
		#ifdef DEBUG
		this->check_bounds(i,j);
		return data_.at(index(i,j,dir));
		#else
		return data_.data()[index(i,j,dir)];
		#endif
	}

	real_t & operator()(len_t i, len_t j, int dir)
	{
		#ifdef DEBUG
		this->check_bounds(i,j);
		return data_.at(index(i,j,dir));
		#else
		return data_.data()[index(i,j,dir)];
		#endif
	}
};


}}}

#endif
