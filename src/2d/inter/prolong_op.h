#ifndef BOXMG_2D_INTER_PROLONG_OP_H
#define BOXMG_2D_INTER_PROLONG_OP_H

#include <utility>
#include <iostream>

#include "boxmg-common.h"

#include "core/stencil_op.h"
#include "core/grid_func.h"


namespace boxmg { namespace bmg2d { namespace inter {

using iadd_pack = std::tuple<const ProlongOp&, const GridFunc &, const GridFunc&>;
class ProlongOp : public StencilOp
{

public:
	ProlongOp() {};
	ProlongOp(len_t nx, len_t ny);
	ProlongOp & operator=(ProlongOp&& P)
	{
		gs = std::move(P.gs);
		return *this;
	}
	ProlongOp(ProlongOp&& P)
	{
		*this = std::move(P);
	}
	friend std::ostream & operator<< (std::ostream &os, const ProlongOp & P);
	friend iadd_pack operator*(const ProlongOp&, const GridFunc&);
	StencilOp * fine_op;
	GridFunc *residual;
};

}}}

#endif
