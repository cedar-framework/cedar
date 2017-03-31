#ifndef CEDAR_2D_CORE_RELAX_STENCIL_H
#define CEDAR_2D_CORE_RELAX_STENCIL_H

#include <cedar/2d/types.h>
#include <cedar/array.h>
#include <cedar/grid_quantity.h>

namespace cedar { namespace cdr2 {


class relax_stencil : public array<len_t, real_t,3>, public grid_quantity<len_t, 2>
{
public:
	relax_stencil() {};
	relax_stencil(len_t nx, len_t ny, unsigned int nghosts=1);

	using array<len_t,real_t,3>::index;
	using array<len_t,real_t,3>::operator();

	/* len_t index(len_t i, len_t j, int dir) const */
	/* { */
	/* 	return dir*len_[0]*len_[1] + i*stride_[0] + j * stride_[1]; */
	/* } */

	/* real_t operator()(len_t i, len_t j, int dir) const */
	/* { */
	/* 	#ifdef DEBUG */
	/* 	this->check_bounds(i,j); */
	/* 	return data_.at(index(i,j,dir)); */
	/* 	#else */
	/* 	return data_.data()[index(i,j,dir)]; */
	/* 	#endif */
	/* } */

	/* real_t & operator()(len_t i, len_t j, int dir) */
	/* { */
	/* 	#ifdef DEBUG */
	/* 	this->check_bounds(i,j); */
	/* 	return data_.at(index(i,j,dir)); */
	/* 	#else */
	/* 	return data_.data()[index(i,j,dir)]; */
	/* 	#endif */
	/* } */
};


}}

#endif
