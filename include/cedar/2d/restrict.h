#ifndef CEDAR_2D_RESTRICTION_H
#define CEDAR_2D_RESTRICTION_H

#include <cedar/kernels/restrict.h>
#include <cedar/2d/types.h>

namespace cedar { namespace cdr2 {

class restrict_f90 : public kernels::restriction<stypes>
{
	void run(const restrict_op & R,
	         const grid_func & x,
	         grid_func & y) override;
};

}}
#endif
