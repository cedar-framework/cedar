#ifndef CEDAR_2D_SOLVE_CG_H
#define CEDAR_2D_SOLVE_CG_H

#include <cedar/2d/types.h>
#include <cedar/kernels/solve_cg.h>


extern "C" {
	using namespace cedar;
	void BMG_get_bc(int, int*);
}


namespace cedar { namespace cdr2 {

class solve_cg_f90 : public kernels::solve_cg<stypes>
{
public:

	solve_cg_f90(kmode kernmode);

	void setup(const stencil_op<five_pt> & so,
	           grid_func & ABD) override
	{
		this->setup_impl(so, ABD);
	}
	void setup(const stencil_op<nine_pt> & so,
	           grid_func & ABD) override
	{
		this->setup_impl(so, ABD);
	}

	template<class sten>
	void setup_impl(const stencil_op<sten> & so,
	                grid_func & ABD)
	{
		len_t nx, ny;
		int nstencil;
		len_t nabd1, nabd2;
		int ibc;

		auto & sod = const_cast<stencil_op<sten>&>(so);

		nx = so.len(0);
		ny = so.len(1);

		nstencil = stencil_ndirs<sten>::value;

		nabd1 = ABD.len(0);
		nabd2 = ABD.len(1);

		BMG_get_bc(params->per_mask(), &ibc);

		fcall_setup(sod.data(), &nx, &ny, &nstencil,
		            ABD.data(), &nabd1, &nabd2, &ibc);
	}

	void run(grid_func & x,
	         const grid_func & b,
	         const grid_func & ABD,
	         real_t * bbd) override;
protected:
	std::function<void(real_t*, real_t*, len_t, len_t, real_t*, real_t*, len_t, len_t, int)> fcall;
	std::function<void(real_t*, len_t*, len_t*, int*, real_t*, len_t*,len_t*,int*)> fcall_setup;
};

}}

#endif

