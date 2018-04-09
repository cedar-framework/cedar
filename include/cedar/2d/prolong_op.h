#ifndef CEDAR_2D_INTER_PROLONG_OP_H
#define CEDAR_2D_INTER_PROLONG_OP_H

#include <utility>
#include <iostream>

#include <cedar/types.h>
#include <cedar/2d/stencil_op.h>
#include <cedar/2d/grid_func.h>


namespace cedar { namespace cdr2 {

enum class inter_dir {L=0,R=1,A=2,B=3,SW=4,NW=5,NE=6,SE=7,ndirs};

class prolong_op : public stencil_op<inter_dir>
{

public:
	prolong_op() {};
	prolong_op(len_t nx, len_t ny);
	friend std::ostream & operator<< (std::ostream &os, const prolong_op & P);
	stencil_op<five_pt> * fine_op_five;
	stencil_op<nine_pt> * fine_op_nine;
	grid_func *residual;
	bool fine_is_five;
};

}}

#endif
