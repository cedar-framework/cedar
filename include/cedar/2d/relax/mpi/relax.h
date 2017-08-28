#ifndef CEDAR_2D_KERNEL_RELAX_MPI_H
#define CEDAR_2D_KERNEL_RELAX_MPI_H

#include "cedar/2d/ftn/BMG_parameters_c.h"
#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"
#include <cedar/kernel_params.h>
#include "cedar/2d/mpi/halo.h"
#include <cedar/halo_exchanger.h>
#include "cedar/cycle/types.h"
#include "cedar/2d/stencil_op.h"
#include "cedar/2d/relax_stencil.h"
#include "cedar/2d/grid_func.h"
#include "cedar/2d/mpi/grid_func.h"
#include "cedar/2d/mpi/stencil_op.h"

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_relax_GS(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                              len_t II, len_t JJ, int kf, int ifd, int nstncl,
	                              int irelax_sym,
	                              int updown, len_t iGs, len_t jGs, void *halof);
	void MPI_BMG2_SymStd_relax_lines_x(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                                   real_t *B, len_t II, len_t JJ, len_t iGs, len_t jGs,
	                                   int nog, int nstencil, int irelax_sym, int updown,
	                                   len_t *datadist, len_t *iwork, len_t nmsgi, int *pMSG,
	                                   real_t *rwork, len_t nmsgr, int mpicomm,
	                                   int xlinecomm, int ylinecomm, void *halof);
	void MPI_BMG2_SymStd_relax_lines_y(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                                   real_t *B, len_t II, len_t JJ, len_t iGs, len_t jGs,
	                                   int nog, int nstencil, int irelax_sym, int updown,
	                                   len_t *datadist, len_t *iwork, len_t nmsgi, int *pMSG,
	                                   real_t *rwork, len_t nmsgr, int mpicomm,
	                                   int xlinecomm, int ylinecomm, void *halof);
}


namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr2::mpi;
	template<class sten>
	void mpi_relax_rbgs_point(const kernel_params & params,
	                          halo_exchanger *halof,
	                          const mpi::stencil_op<sten> & so,
	                          mpi::grid_func & x,
	                          const mpi::grid_func & b,
	                          const relax_stencil & sor,
	                          cycle::Dir cycle_dir)
	{
		using namespace cedar::cdr2;
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

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		// ibc = BMG_BCs_definite;

		MPI_BMG2_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		                         so.len(0), so.len(1), kf, ifd, nstencil, BMG_RELAX_SYM,
		                         updown, topo.is(0), topo.is(1), halof);
	}

	template<class sten>
	void mpi_relax_lines_x(const kernel_params & params,
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
		MsgCtx *ctx = (MsgCtx*) halof->context_ptr();
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		mpi::grid_func & bd = const_cast<mpi::grid_func&>(b);

		k = topo.level()+1;
		kf = topo.nlevel();

		nstencil = stencil_ndirs<sten>::value;

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		// ibc = BMG_BCs_definite;
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		MPI_Fint xlinecomm = MPI_Comm_c2f(ctx->xlinecomm);
		MPI_Fint ylinecomm = MPI_Comm_c2f(ctx->ylinecomm);

		cedar::len_t * xdatadist = &ctx->msg_geom.data()[ctx->pLS(ipL_LS_XDataDist,k-1)-1];

		MPI_BMG2_SymStd_relax_lines_x(k, sod.data(), bd.data(), x.data(), sord.data(), res.data(),
		                              so.len(0), so.len(1), topo.is(0), topo.is(1),
		                              kf, nstencil, BMG_RELAX_SYM, updown,
		                              xdatadist,
		                              ctx->msg_geom.data(), ctx->msg_geom.size(),
		                              ctx->pMSG.data(), ctx->msg_buffer.data(),
		                              ctx->msg_buffer.size(), fcomm,
		                              xlinecomm, ylinecomm, halof);
	}

	template<class sten>
	void mpi_relax_lines_y(const kernel_params & params,
	                       halo_exchanger *halof,
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
		MsgCtx *ctx = (MsgCtx*) halof->context_ptr();
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		mpi::grid_func & bd = const_cast<mpi::grid_func&>(b);

		k = topo.level()+1;
		kf = topo.nlevel();

		nstencil = stencil_ndirs<sten>::value;

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		// ibc = BMG_BCs_definite;
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		MPI_Fint xlinecomm = MPI_Comm_c2f(ctx->xlinecomm);
		MPI_Fint ylinecomm = MPI_Comm_c2f(ctx->ylinecomm);

		cedar::len_t * ydatadist = &ctx->msg_geom.data()[ctx->pLS(ipL_LS_YDataDist,k-1)-1];

		MPI_BMG2_SymStd_relax_lines_y(k, sod.data(), bd.data(), x.data(), sord.data(), res.data(),
		                              so.len(0), so.len(1), topo.is(0), topo.is(1),
		                              kf, nstencil, BMG_RELAX_SYM, updown,
		                              ydatadist,
		                              ctx->msg_geom.data(), ctx->msg_geom.size(),
		                              ctx->pMSG.data(), ctx->msg_buffer.data(),
		                              ctx->msg_buffer.size(), fcomm,
		                              xlinecomm, ylinecomm, halof);
	}
}

}}}

#endif
