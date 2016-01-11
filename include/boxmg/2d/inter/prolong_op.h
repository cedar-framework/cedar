#ifndef BOXMG_2D_INTER_PROLONG_OP_H
#define BOXMG_2D_INTER_PROLONG_OP_H

#include <utility>
#include <iostream>

#include <boxmg/types.h>
#include <boxmg/2d/stencil_op.h>
#include <boxmg/2d/grid_func.h>


namespace boxmg { namespace bmg2d { namespace inter {

using iadd_pack = std::tuple<const prolong_op&, const grid_func &, const grid_func&>;
class prolong_op : public stencil_op
{

public:
	prolong_op() {};
	prolong_op(len_t nx, len_t ny);
	prolong_op & operator=(prolong_op&& P)
	{
		gs = std::move(P.gs);
		return *this;
	}
	prolong_op(prolong_op&& P)
	{
		*this = std::move(P);
	}
	friend std::ostream & operator<< (std::ostream &os, const prolong_op & P);
	friend iadd_pack operator*(const prolong_op&, const grid_func&);
	stencil_op * fine_op;
	grid_func *residual;
};

}}}

#endif
