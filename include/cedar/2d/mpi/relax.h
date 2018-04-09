#ifndef CEDAR_2D_RELAX_MPI_H
#define CEDAR_2D_RELAX_MPI_H

#include <cedar/2d/ftn/mpi/BMG_parameters_c.h>
#include <cedar/2d/ftn/mpi/BMG_workspace_c.h>

#include <cedar/2d/mpi/types.h>
#include <cedar/kernels/point_relax.h>
#include <cedar/kernels/line_relax.h>

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

		MPI_BMG2_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		                         so.len(0), so.len(1), kf, ifd, nstencil, BMG_RELAX_SYM,
		                         updown, topo.is(0), topo.is(1), halof);
	}
};


template<relax_dir rdir>
class lines : public kernels::line_relax<mpi::stypes, rdir>
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
		int nstencil = stencil_ndirs<sten>::value;
		auto & sod = const_cast<stencil_op<sten>&>(so);

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
		MPI_Fint xlinecomm = MPI_Comm_c2f(this->halof->linecomm(0));
		MPI_Fint ylinecomm = MPI_Comm_c2f(this->halof->linecomm(1));

		cedar::len_t * xdatadist = this->halof->datadist(0,k);
		cedar::len_t * ydatadist = this->halof->datadist(1,k);

		if (rdir == relax_dir::x) {
			MPI_BMG2_SymStd_relax_lines_x(k, sod.data(), bd.data(), x.data(), sord.data(), res.data(),
			                              so.len(0), so.len(1), topo.is(0), topo.is(1),
			                              kf, nstencil, BMG_RELAX_SYM, updown,
			                              xdatadist,
			                              this->halof->linebuf().data(),
			                              this->halof->linebuf().size(), fcomm,
			                              xlinecomm, ylinecomm, this->halof);
		} else {
			MPI_BMG2_SymStd_relax_lines_y(k, sod.data(), bd.data(), x.data(), sord.data(), res.data(),
			                              so.len(0), so.len(1), topo.is(0), topo.is(1),
			                              kf, nstencil, BMG_RELAX_SYM, updown,
			                              ydatadist,
			                              this->halof->linebuf().data(),
			                              this->halof->linebuf().size(), fcomm,
			                              xlinecomm, ylinecomm, this->halof);
		}
	}
};

}}}
#endif
