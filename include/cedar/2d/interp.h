#ifndef CEDAR_2D_INTERP_H
#define CEDAR_2D_INTERP_H


#include <cedar/2d/types.h>
#include <cedar/2d/prolong_op.h>
#include <cedar/2d/grid_func.h>
#include <cedar/kernels/interp_add.h>
#include <cedar/kernels/setup_interp.h>

namespace cedar { namespace cdr2 {

class interp_f90 : public kernels::interp_add<cdr2::stypes>
{
public:
	using prolong_op = cedar::cdr2::prolong_op;
	using grid_func = cedar::cdr2::grid_func;
	interp_f90(kmode kernmode);
	void run(const prolong_op & P,
	         const grid_func & coarse,
	         const grid_func & residual,
	         grid_func & fine) override;
protected:
	std::function<void(real_t*, real_t*, real_t*,real_t*,real_t*,
	                   len_t, len_t, len_t, len_t, int, int)> fcall;
};


class setup_interp_f90 : public kernels::setup_interp<cdr2::stypes>
{
	void run(const stencil_op<five_pt> & fop,
	         const stencil_op<nine_pt> & cop,
	         prolong_op & P) override;
	void run(const stencil_op<nine_pt> & fop,
	         const stencil_op<nine_pt> & cop,
	         prolong_op & P) override;
};

}}

#endif
