#ifndef CEDAR_2D_RELAX_MPI_H
#define CEDAR_2D_RELAX_MPI_H

#include <cedar/2d/ftn/mpi/BMG_parameters_c.h>
#include <cedar/2d/ftn/mpi/BMG_workspace_c.h>
#include <cedar/device.h>

#include <cedar/2d/mpi/types.h>
#include <cedar/kernels/point_relax.h>
#include <cedar/kernels/line_relax.h>
#include <cedar/2d/mpi/kernel_manager.h>

using real_t = cedar::real_t;
using len_t = cedar::len_t;
#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_recip.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_relax_GS.f90.hpp>

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_SETUP_recip(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
	void MPI_BMG2_SymStd_SETUP_lines_x(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
	void MPI_BMG2_SymStd_SETUP_lines_y(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
	void MPI_BMG2_SymStd_relax_GS(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                              len_t II, len_t JJ, int kf, int ifd, int nstncl,
	                              int irelax_sym,
	                              int updown, len_t iGs, len_t jGs, void *halof);
	void MPI_BMG2_SymStd_relax_lines_x(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                                   real_t *B, len_t II, len_t JJ, len_t iGs, len_t jGs,
	                                   int nog, int nstencil, int irelax_sym, int updown,
	                                   len_t *datadist,
	                                   real_t *rwork, len_t nmsgr, int mpicomm,
	                                   int xlinecomm, int ylinecomm, void *halof);
	void MPI_BMG2_SymStd_relax_lines_y(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                                   real_t *B, len_t II, len_t JJ, len_t iGs, len_t jGs,
	                                   int nog, int nstencil, int irelax_sym, int updown,
	                                   len_t *datadist,
	                                   real_t *rwork, len_t nmsgr, int mpicomm,
	                                   int xlinecomm, int ylinecomm, void *halof);
}


namespace cedar { namespace cdr2 { namespace mpi {

template <typename device=cedar::cpu>
class rbgs : public kernels::point_relax<mpi::stypes>
{
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

                sod.template ensure<device>();
                sor.template ensure<device>();

                if (device::is_gpu()) {
                    MPI_BMG2_SymStd_SETUP_recip<ftl::device::GPU>(sod, sor, so.len(0), so.len(1), nstencil);
                } else {
                    MPI_BMG2_SymStd_SETUP_recip(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
                    sor.template mark_dirty<device>();
                }
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

                sod.template ensure<device>();
                x.template ensure<device>();
                bd.template ensure<device>();
                sord.template ensure<device>();

                void* halof = this->services->fortran_handle<halo_exchange>();

                if (device::is_gpu()) {
                    MPI_BMG2_SymStd_relax_GS<cedar::gpu>(
                        k, sod, bd, x, sord, so.len(0), so.len(1), kf, ifd, nstencil,
                        BMG_RELAX_SYM, updown, topo.is(0), topo.is(1), halof);
                } else {
                    MPI_BMG2_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
                                             so.len(0), so.len(1), kf, ifd, nstencil, BMG_RELAX_SYM,
                                             updown, topo.is(0), topo.is(1), halof);
                    x.template mark_dirty<device>();
                }
	}
};


template<relax_dir rdir>
class lines : public kernels::line_relax<mpi::stypes, rdir>
{
public:
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
		int nstencil = stencil_ndirs<sten>::value;
		auto & sod = const_cast<stencil_op<sten>&>(so);

		auto & topo = so.grid();
		auto & mp = this->services->template get<message_passing>();
		mp.comm_split(topo.comm, topo.coord(1), topo.coord(0), &xlinecomm);
		mp.comm_split(topo.comm, topo.coord(0), topo.coord(1), &ylinecomm);

		if (rdir == relax_dir::x)
			MPI_BMG2_SymStd_SETUP_lines_x(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
		else
			MPI_BMG2_SymStd_SETUP_lines_y(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
	}

	template<class sten>
	void run_impl(const stencil_op<sten> & so,
	              grid_func & x,
	              const grid_func & b,
	              const relax_stencil & sor,
	              grid_func & res,
	              cycle::Dir cdir)
	{
		using namespace cedar::cdr2;
		int k, kf;
		int updown, nstencil;

		auto & sod = const_cast<mpi::stencil_op<sten>&>(so);
		grid_topo & topo = sod.grid();
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		mpi::grid_func & bd = const_cast<mpi::grid_func&>(b);

		k = topo.level()+1;
		kf = topo.nlevel();

		nstencil = stencil_ndirs<sten>::value;

		if (cdir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		// ibc = BMG_BCs_definite;
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		MPI_Fint fxlinecomm = MPI_Comm_c2f(this->xlinecomm);
		MPI_Fint fylinecomm = MPI_Comm_c2f(this->ylinecomm);

		auto & halo_service = this->services->template get<halo_exchange>();
		cedar::len_t * xdatadist = halo_service.datadist(0,k);
		cedar::len_t * ydatadist = halo_service.datadist(1,k);

		if (rdir == relax_dir::x) {
			MPI_BMG2_SymStd_relax_lines_x(k, sod.data(), bd.data(), x.data(), sord.data(), res.data(),
			                              so.len(0), so.len(1), topo.is(0), topo.is(1),
			                              kf, nstencil, BMG_RELAX_SYM, updown,
			                              xdatadist,
			                              halo_service.linebuf().data(),
			                              halo_service.linebuf().size(), fcomm,
			                              fxlinecomm, fylinecomm,
			                              this->services->template fortran_handle<halo_exchange>());
		} else {
			MPI_BMG2_SymStd_relax_lines_y(k, sod.data(), bd.data(), x.data(), sord.data(), res.data(),
			                              so.len(0), so.len(1), topo.is(0), topo.is(1),
			                              kf, nstencil, BMG_RELAX_SYM, updown,
			                              ydatadist,
			                              halo_service.linebuf().data(),
			                              halo_service.linebuf().size(), fcomm,
			                              fxlinecomm, fylinecomm,
			                              this->services->template fortran_handle<halo_exchange>());
		}
	}


protected:
	MPI_Comm xlinecomm;
	MPI_Comm ylinecomm;
};

}}}
#endif
