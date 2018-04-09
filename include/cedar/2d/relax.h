#ifndef CEDAR_2D_RELAX_SER_H
#define CEDAR_2D_RELAX_SER_H

#include <cedar/types.h>
#include <cedar/2d/types.h>
#include <cedar/kernels/point_relax.h>
#include <cedar/kernels/line_relax.h>
#include <cedar/2d/ftn/BMG_parameters_c.h>


extern "C" {
	using namespace cedar;
	void BMG2_SymStd_SETUP_recip(real_t *so, real_t *sor, len_t nx, len_t ny, int nstncl, int nsor_v);
	void BMG2_SymStd_SETUP_lines_x(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl, int JPN);
	void BMG2_SymStd_SETUP_lines_y(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl, int JPN);
	void BMG2_SymStd_relax_GS(int, real_t*, real_t*, real_t*, real_t*, len_t, len_t,
	                          int, int, int, int, int, int, int);
	void BMG2_SymStd_relax_lines_x(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                               real_t *B, len_t II, len_t JJ, int kf, int ifd,
	                               int nstencil, int irelax_sym, int updown, int jpn);
	void BMG2_SymStd_relax_lines_y(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                               real_t *B, len_t II, len_t JJ, int kf, int ifd,
	                               int nstencil, int irelax_sym, int updown, int jpn);
	void BMG_get_bc(int, int*);
}


namespace cedar { namespace cdr2 {

class rbgs : public kernels::point_relax<stypes>
{
	virtual void setup(const stencil_op<five_pt> & so,
	                   relax_stencil & sor) override
	{
		this->setup_impl(so, sor);
	}
	virtual void setup(const stencil_op<nine_pt> & so,
	                   relax_stencil & sor) override
	{
		this->setup_impl(so, sor);
	}
	virtual void run(const stencil_op<five_pt> & so,
	                 grid_func & x,
	                 const grid_func & b,
	                 const relax_stencil & sor,
	                 cycle::Dir cdir) override
	{
		this->run_impl(so, x, b, sor, cdir);
	}
	virtual void run(const stencil_op<nine_pt> & so,
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
		int nstencil = stencil_ndirs<sten>::value;
		int nsorv = 2;

		auto & sod = const_cast<stencil_op<sten>&>(so);

		BMG2_SymStd_SETUP_recip(sod.data(), sor.data(), so.len(0), so.len(1), nstencil, nsorv);
	}


	template<class sten>
	void run_impl(const stencil_op<sten> & so,
	              grid_func & x,
	              const grid_func & b,
	              const relax_stencil & sor,
	              cycle::Dir cdir)
	{
		int k, kf, ifd;
		int updown, nsorv, ibc, nstencil;

		auto & sod = const_cast<stencil_op<sten>&>(so);
		auto & bd = const_cast<grid_func&>(b);
		auto & sord = const_cast<relax_stencil&>(sor);
		k = kf = 1;
		if (std::is_same<sten, five_pt>::value)
			ifd = 1;
		else
			ifd = 0;

		nstencil = stencil_ndirs<sten>::value;
		nsorv = 2;

		if (cdir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		BMG_get_bc(params->per_mask(), &ibc);

		BMG2_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		                     so.len(0), so.len(1), kf, ifd, nstencil, nsorv,
		                     BMG_RELAX_SYM, updown, ibc);
	}
};


template<relax_dir rdir>
class lines : public kernels::line_relax<stypes, rdir>
{
	virtual void setup(const stencil_op<five_pt> & so,
	                   relax_stencil & sor) override
	{
		this->setup_impl(so, sor);
	}
	virtual void setup(const stencil_op<nine_pt> & so,
	                   relax_stencil & sor) override
	{
		this->setup_impl(so, sor);
	}
	virtual void run(const stencil_op<five_pt> & so,
	                 grid_func & x,
	                 const grid_func & b,
	                 const relax_stencil & sor,
	                 grid_func & res,
	                 cycle::Dir cdir) override
	{
		this->run_impl(so, x, b, sor, res, cdir);
	}
	virtual void run(const stencil_op<nine_pt> & so,
	                 grid_func & x,
	                 const grid_func & b,
	                 const relax_stencil & sor,
	                 grid_func & res,
	                 cycle::Dir cdir) override
	{
		this->run_impl(so, x, b, sor, res, cdir);
	}


	template<class sten>
	void setup_impl(const stencil_op<sten> & so,
	                relax_stencil & sor)
	{
		int jpn;
		int nstencil = stencil_ndirs<sten>::value;

		auto & sod = const_cast<stencil_op<sten>&>(so);

		BMG_get_bc(this->params->per_mask(), &jpn);

		std::function<void(real_t*, real_t*, len_t, len_t, int, int)> relax;

		if (rdir == relax_dir::x)
			BMG2_SymStd_SETUP_lines_x(sod.data(), sor.data(), so.len(0), so.len(1), nstencil, jpn);
		else
			BMG2_SymStd_SETUP_lines_y(sod.data(), sor.data(), so.len(0), so.len(1), nstencil, jpn);
	}


	template<class sten>
	void run_impl(const stencil_op<sten> & so,
	              grid_func & x,
	              const grid_func & b,
	              const relax_stencil & sor,
	              grid_func & res,
	              cycle::Dir cdir)
	{
		int k, kf, ifd;
		int updown, ibc, nstencil;

		auto & sod = const_cast<stencil_op<sten>&>(so);
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		grid_func & bd = const_cast<grid_func&>(b);
		k = kf = 1;
		nstencil = stencil_ndirs<sten>::value;
		if (std::is_same<sten, five_pt>::value)
			ifd = 1;
		else
			ifd = 0;

		if (cdir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		BMG_get_bc(this->params->per_mask(), &ibc);

		if (rdir == relax_dir::x) {
			BMG2_SymStd_relax_lines_x(k, sod.data(), bd.data(), x.data(), sord.data(), res.data(),
			                          so.len(0), so.len(1), kf, ifd, nstencil,
			                          BMG_RELAX_SYM, updown, ibc);
		} else {
			BMG2_SymStd_relax_lines_y(k, sod.data(), bd.data(), x.data(), sord.data(), res.data(),
			                          so.len(0), so.len(1), kf, ifd, nstencil,
			                          BMG_RELAX_SYM, updown, ibc);
		}
	}

};

}}

#endif
