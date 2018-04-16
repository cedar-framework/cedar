#ifndef CEDAR_3D_MPI_RELAX_H
#define CEDAR_3D_MPI_RELAX_H

#include <cedar/2d/ftn/mpi/BMG_parameters_c.h>
#include <cedar/3d/mpi/types.h>
#include <cedar/kernels/point_relax.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG3_SymStd_SETUP_recip(real_t *so, real_t *sor,
	                                 len_t nx, len_t ny, len_t nz,
	                                 int nstencil, int nsorv);
	void MPI_BMG3_SymStd_relax_GS(int kg, real_t *so, real_t *qf, real_t *q, real_t *sor,
	                              len_t nlx, len_t nly, len_t nlz, len_t ngx, len_t ngy, len_t ngz,
	                              int nog, int ifd, int nstencil, int nsorv,
	                              int irelax_sym, int updown,
	                              len_t igs, len_t jgs, len_t kgs, void *halof);
}

namespace cedar { namespace cdr3 { namespace mpi {

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
		int nstencil, nsorv;

		auto & sod = const_cast<stencil_op<sten>&>(so);

		nstencil = stencil_ndirs<sten>::value;

		nsorv = 2;

		MPI_BMG3_SymStd_SETUP_recip(sod.data(),
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
		int k, ifd;
		int updown, nstencil;

		auto & sod = const_cast<mpi::stencil_op<sten>&>(so);
		grid_topo & topo = sod.grid();
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		mpi::grid_func & bd = const_cast<mpi::grid_func&>(b);

		k = topo.level()+1;
		nstencil = stencil_ndirs<sten>::value;

		if (std::is_same<sten, seven_pt>::value)
			ifd = 1;
		else
			ifd = 0;

		if (cdir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		int nsorv = 2;

		// ibc = BMG_BCs_definite;
		MPI_BMG3_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		                         so.len(0), so.len(1), so.len(2),
		                         topo.nglobal(0), topo.nglobal(1), topo.nglobal(2),
		                         topo.nlevel(), ifd, nstencil, nsorv,
		                         BMG_RELAX_SYM, updown,
		                         topo.is(0), topo.is(1), topo.is(2),
		                         halof);
	}
};
}}}

#endif

