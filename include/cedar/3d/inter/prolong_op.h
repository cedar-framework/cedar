#ifndef CEDAR_3D_INTER_PROLONG_OP_H
#define CEDAR_3D_INTER_PROLONG_OP_H

#include <utility>
#include <iostream>

#include <cedar/types.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/grid_func.h>

namespace cedar { namespace cdr3 { namespace inter {

using iadd_pack = std::tuple<const prolong_op&, const grid_func&, const grid_func&>;
class prolong_op : public stencil_op
{
public:
	prolong_op() {}
	prolong_op(len_t nx, len_t ny, len_t nz);
	prolong_op & operator=(prolong_op&& P)
	{
		gs = std::move(P.gs);
		return *this;
	}
	prolong_op(prolong_op&& P)
	{
		*this = std::move(P);
	}
	friend std::ostream &operator<< (std::ostream & os, const prolong_op &P);
	friend iadd_pack operator*(const prolong_op&, const grid_func&);
	stencil_op *fine_op;
	grid_func *residual;
};

}}}


#endif
