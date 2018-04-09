#ifndef CEDAR_3D_INTER_PROLONG_OP_H
#define CEDAR_3D_INTER_PROLONG_OP_H

#include <utility>
#include <iostream>

#include <cedar/types.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/grid_func.h>

namespace cedar { namespace cdr3 {

enum class inter_dir {XYL=0, XYR=1, XYA=2, XYB=3, XZA=4, XZB=5,
                      XYNE=6, XYSE=7, XYSW=8, XYNW=9, XZSW=10,
                      XZNW=11,
                      XZNE=12, XZSE=13, YZSW=14, YZNW=15, YZNE=16,
                      YZSE=17, BSW=18, BNW=19, BNE=20, BSE=21,
                      TSW=22, TNW=23, TNE=24, TSE=25, ndirs};

class prolong_op : public stencil_op<inter_dir>
{
public:
	prolong_op() {}
	prolong_op(len_t nx, len_t ny, len_t nz);
	friend std::ostream &operator<< (std::ostream & os, const prolong_op &P);
	stencil_op<seven_pt> * fine_op_seven;
	stencil_op<xxvii_pt> * fine_op_xxvii;
	grid_func * residual;
	bool fine_is_seven;
};

}}


#endif
