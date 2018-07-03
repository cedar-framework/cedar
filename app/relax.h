#ifndef CEDAR_APP_RELAX_H
#define CEDAR_APP_RELAX_H

#include <cedar/2d/ftn/mpi/BMG_parameters_c.h>
#include <cedar/2d/ftn/mpi/BMG_workspace_c.h>
#include <cedar/2d/mpi/types.h>
#include <cedar/kernels/point_relax.h>

extern "C" {
	using namespace cedar;

	void MPI_BMG2_SymStd_relax_GS_plane(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                                    len_t II, len_t JJ, int nplanes, int kf, int ifd, int nstncl,
	                                    int irelax_sym,
	                                    int updown, len_t iGs, len_t jGs, void *halof);
	void MPI_BMG2_SymStd_SETUP_recip(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
}


namespace cedar { namespace cdr2 { namespace mpi {

class rbgs_plane : public kernels::point_relax<mpi::stypes>
{
public:
	rbgs_plane(int nplanes) : nplanes(nplanes) {}
	void setup(const stencil_op<five_pt> & so,
	           relax_stencil & sor) override
	{
		this->setup_impl(so, sor);
	}
	void setup(const stencil_op<nine_pt> & so,
	           relax_stencil & sor) override
	{
		this->setup_impl(so, sor);
	}
	void run(const stencil_op<five_pt> & so,
	         grid_func & x,
	         const grid_func & b,
	         const relax_stencil & sor,
	         cycle::Dir cdir) override
	{
		this->run_impl(so, x, b, sor, cdir);
	}
	void run(const stencil_op<nine_pt> & so,
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
		auto & sod = const_cast<stencil_op<sten>&>(so);
		MPI_BMG2_SymStd_SETUP_recip(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
	}

	template<class sten>
	void run_impl(const stencil_op<sten> & so,
	              grid_func & x,
	              const grid_func & b,
	              const relax_stencil & sor,
	              cycle::Dir cdir)
	{
		int k, kf, ifd;
		int updown, nstencil;

		auto & sod = const_cast<mpi::stencil_op<sten>&>(so);
		grid_topo & topo = sod.grid();
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		mpi::grid_func & bd = const_cast<mpi::grid_func&>(b);

		k = topo.level()+1;
		kf = topo.nlevel();

		if (std::is_same<five_pt, sten>::value)
			ifd = 1;
		else
			ifd = 0;
		nstencil = stencil_ndirs<sten>::value;

		if (cdir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		// ibc = BMG_BCs_definite;

		MPI_BMG2_SymStd_relax_GS_plane(k, sod.data(), bd.data(), x.data(), sord.data(),
		                               so.len(0), so.len(1), this->nplanes, kf, ifd, nstencil, BMG_RELAX_SYM,
		                               updown, topo.is(0), topo.is(1),
		                               services->fortran_handle<halo_exchange>());
	}
protected:
	int nplanes;
};

}}}

#endif
