#ifndef CEDAR_2D_KERNEL_RELAX_MPI_ML_H
#define CEDAR_2D_KERNEL_RELAX_MPI_ML_H

#include "cedar/2d/ftn/BMG_parameters_c.h"
#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"
#include <cedar/2d/mpi/msg_exchanger.h>
#include <cedar/2d/mpi/tausch_exchanger.h>


extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_SETUP_trisolve(int nproc, len_t nlines, int min_gsz, int nog, int *nspace,
	                               int *nol);
	void BMG2_SymStd_SETUP_ADD_SOR_PTRS(int nproc, len_t nlines, int min_gsz, int nol, int *poffset,
	                                    int *psor_lev);
	void MPI_BMG2_SymStd_SETUP_comm(MPI_Fint *comm, int nol, int rank, int min_gsz);
	void MPI_BMG2_SymStd_SETUP_lines_x(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
	void MPI_BMG2_SymStd_SETUP_lines_y(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
	void MPI_BMG2_SymStd_relax_lines_x_ml(int k, real_t *so, real_t *qf, real_t *q,
	                                      real_t *sor, real_t *res,
	                                      len_t II, len_t JJ, len_t igs, len_t jgs,
	                                      int nog, int nstencil, int irelax_sym, int updown,
	                                      len_t *datadist, real_t *rwork, len_t nmsgr,
	                                      MPI_Fint mpicomm, MPI_Fint *xcomm, int nolx,
	                                      int *tdsx_sor_ptrs, len_t nsor, real_t *tdg, void *halof);
}

namespace cedar { namespace cdr2 { namespace mpi {


enum class relax_dir { x, y };

template<class halo_exchanger>
	len_t * get_datadist(halo_exchanger * halof, int k, relax_dir dir);

template<>
	len_t * get_datadist(mpi::msg_exchanger * halof, int k, relax_dir dir)
{
	len_t * ret;
	kernel::impls::MsgCtx * ctx = (kernel::impls::MsgCtx*) halof->context_ptr();

	if (dir == relax_dir::x)
		ret = &ctx->msg_geom.data()[ctx->pLS(ipL_LS_XDataDist,k-1)-1];
	else
		ret = &ctx->msg_geom.data()[ctx->pLS(ipL_LS_YDataDist,k-1)-1];

	return ret;
}

template<>
	len_t * get_datadist(mpi::tausch_exchanger * halof, int k, relax_dir dir)
{
	len_t * ret;
	mpi::line_pkg & line_info = halof->line_data;

	if (dir == relax_dir::x)
		ret = line_info.datadist[0].data();
	else
		ret = line_info.datadist[1].data();

	return ret;
}


class ml_line_relaxer
{
public:
	ml_line_relaxer(relax_dir dir) : dir(dir) {}

	void init_ndim(int nproc, int nproc_other, int coord, int coord_other, len_t nlines, len_t nlines_other,
	               int min_gsz, int nog, MPI_Comm comm)
	{
			int workspace_size, nspace;
			MPI_BMG2_SymStd_SETUP_trisolve(nproc, nlines, min_gsz, nog,
			                               &nspace, &this->nlevels);
			comms.init(2, this->nlevels - 1);
			MPI_Comm gcomm;
			MPI_Comm_split(comm, coord_other, coord, &gcomm);
			comms(0, this->nlevels - 2) = MPI_Comm_c2f(gcomm);
			MPI_BMG2_SymStd_SETUP_comm(comms.data(), this->nlevels, coord+1, min_gsz);
			sor_ptrs.init(this->nlevels);
			BMG2_SymStd_SETUP_ADD_SOR_PTRS(nproc, nlines, min_gsz, this->nlevels,
			                               &workspace_size, sor_ptrs.data());
			workspace_size += (nlines + 2)*(nlines_other + 2)*4;
			tdg.init(workspace_size);
			int rwork_size = std::max(4 * nlines_other * nproc + 5 * nlines,
			                          4 * nlines * nproc_other + 5 * nlines_other);
			rwork.init(rwork_size);
	}




	template<class sten>
		void setup(const kernel_params & params,
		           const mpi::stencil_op<sten> & so,
		           relax_stencil & sor)
	{
		int nstencil = stencil_ndirs<sten>::value;
		auto & sod = const_cast<mpi::stencil_op<sten>&>(so);

		auto & topo = sod.grid();
		int min_gsz = 3;

		if (dir == relax_dir::x) {
			init_ndim(topo.nproc(0), topo.nproc(1),
			          topo.coord(0), topo.coord(1),
			          topo.nlocal(0) - 2, topo.nlocal(1) - 2,
			          min_gsz, topo.nlevel(), topo.comm);
			MPI_BMG2_SymStd_SETUP_lines_x(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
		} else {
			init_ndim(topo.nproc(1), topo.nproc(0),
			          topo.coord(1), topo.coord(0),
			          topo.nlocal(1) - 2, topo.nlocal(0) - 2,
			          min_gsz, topo.nlevel(), topo.comm);
			MPI_BMG2_SymStd_SETUP_lines_y(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
		}
	}


	template<class sten, class halo_exchanger>
		void solve(const kernel_params & params,
		           halo_exchanger * halof,
		           const mpi::stencil_op<sten> & so,
		           mpi::grid_func & x,
		           const mpi::grid_func & b,
		           const relax_stencil & sor,
		           mpi::grid_func & res,
		           cycle::Dir cycle_dir)
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

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		// ibc = BMG_BCs_definite;
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		len_t * datadist = get_datadist<halo_exchanger>(halof, k, dir);

		std::function<void(int k, real_t *so, real_t *qf, real_t *q,
	                                      real_t *sor, real_t *res,
	                                      len_t II, len_t JJ, len_t igs, len_t jgs,
	                                      int nog, int nstencil, int irelax_sym, int updown,
	                                      len_t *datadist, real_t *rwork, len_t nmsgr,
	                                      MPI_Fint mpicomm, MPI_Fint *xcomm, int nolx,
		                   int *tdsx_sor_ptrs, len_t nsor, real_t *tdg, void *halof)> relax_lines;

		if (dir == relax_dir::x)
			relax_lines = MPI_BMG2_SymStd_relax_lines_x_ml;

		relax_lines(k, sod.data(), bd.data(), x.data(), sord.data(), res.data(),
		            so.len(0), so.len(1), topo.is(0), topo.is(1),
		            kf, nstencil, BMG_RELAX_SYM, updown,
		            datadist, rwork.data(), rwork.len(0),
		            fcomm, this->comms.data(), this->nlevels,
		            sor_ptrs.data(), tdg.len(0), tdg.data(), halof);
	}


protected:
	relax_dir dir;
	int nlevels;
	array<int, 1> sor_ptrs;
	array<real_t, 1> tdg;
	array<real_t, 1> rwork;
	array<MPI_Fint, 2> comms;
};

}}}

#endif
