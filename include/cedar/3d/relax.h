#ifndef CEDAR_3D_RELAX_H
#define CEDAR_3D_RELAX_H

#include <cedar/types.h>
#include <cedar/3d/types.h>
#include <cedar/kernels/point_relax.h>
#include <cedar/2d/ftn/BMG_parameters_c.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SETUP_recip(real_t *so, real_t *sor,
	                             len_t nx, len_t ny, len_t nz,
	                             int nstencl, int nsorv);
	void BMG3_SymStd_relax_GS(int kg, real_t *so, real_t *qf, real_t *q, real_t *sor,
	                          len_t ii, len_t jj, len_t kk, int ifd, int nstncl, int nsorv,
	                          int irelax_sym, int updown, int jpn);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr3 {

class rbgs : public kernels::point_relax<stypes>
{
	void setup(const stencil_op<seven_pt> & so,
	           relax_stencil & sor) override
	{
		this->setup_impl(so, sor);
	}
	void setup(const stencil_op<xxvii_pt> & so,
	           relax_stencil & sor) override
	{
		this->setup_impl(so, sor);
	}
	void run(const stencil_op<seven_pt> & so,
	         grid_func & x,
	         const grid_func & b,
	         const relax_stencil & sor,
	         cycle::Dir cdir) override
	{
		this->run_impl(so, x, b, sor, cdir);
	}
	void run(const stencil_op<xxvii_pt> & so,
	         grid_func & x,
	         const grid_func & b,
	         const relax_stencil & sor,
	         cycle::Dir cdir) override
	{
		this->run_impl(so, x, b, sor, cdir);
	}

	template<class sten>
	void setup_impl(const stencil_op<sten> & so,
	                relax_stencil & sor)
	{
		int nsorv, nstencil;
		auto &sod = const_cast<stencil_op<sten>&>(so);

		nstencil = stencil_ndirs<sten>::value;

		nsorv = 2;

		BMG3_SymStd_SETUP_recip(sod.data(),
		                        sor.data(),
		                        so.len(0), so.len(1), so.len(2),
		                        nstencil, nsorv);
	}

	template<class sten>
	void run_impl(const stencil_op<sten> & so,
	              grid_func & x,
	              const grid_func & b,
	              const relax_stencil & sor,
	              cycle::Dir cdir)
	{
		int k, ifd, nstencil, nsorv, jpn, updown;

		auto & sod = const_cast<stencil_op<sten>&>(so);
		auto & sord = const_cast<relax_stencil&>(sor);
		auto & bd = const_cast<grid_func&>(b);

		nsorv = 2;
		// give these dummy values, ifd will pass the relevant information.
		k = 1;

		nstencil = stencil_ndirs<sten>::value;
		if (std::is_same<sten, seven_pt>::value)
			ifd = 1;
		else
			ifd = 0;

		if (cdir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		BMG_get_bc(params->per_mask(), &jpn);

		BMG3_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		                     so.len(0), so.len(1), so.len(2), ifd, nstencil, nsorv,
		                     BMG_RELAX_SYM, updown, jpn);
	}
};

}}
#endif

