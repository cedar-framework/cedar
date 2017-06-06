#include "cedar/2d/ftn/BMG_parameters_c.h"
#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"
#include "cedar/2d/mpi/stencil_op.h"

#include "cedar/2d/mpi/halo.h"
#include "cedar/2d/relax/relax.h"

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_relax_GS(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                              len_t II, len_t JJ, int kf, int ifd, int nstncl, int irelax_sym,
	                              int updown, len_t iGs, len_t jGs, len_t *iWork, len_t NMSGi,
	                              int *pMSG, real_t *msg_buffer, len_t NMSGr, int MPICOMM);
	void MPI_BMG2_SymStd_relax_lines_x(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                                   real_t *B, len_t II, len_t JJ, len_t iGs, len_t jGs,
	                                   int nog, int nstencil, int irelax_sym, int updown,
	                                   len_t *datadist, len_t *iwork, len_t nmsgi, int *pMSG,
	                                   real_t *rwork, len_t nmsgr, int mpicomm,
	                                   int xlinecomm, int ylinecomm);
	void MPI_BMG2_SymStd_relax_lines_y(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                                   real_t *B, len_t II, len_t JJ, len_t iGs, len_t jGs,
	                                   int nog, int nstencil, int irelax_sym, int updown,
	                                   len_t *datadist, len_t *iwork, len_t nmsgi, int *pMSG,
	                                   real_t *rwork, len_t nmsgr, int mpicomm,
	                                   int xlinecomm, int ylinecomm);
}

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	template<>
	void mpi_relax_rbgs_point(const kernel_params & params,
	                          const mpi::stencil_op<five_pt> & so,
	                          mpi::grid_func & x,
	                          const mpi::grid_func & b,
	                          const relax_stencil & sor,
	                          cycle::Dir cycle_dir)
	{
		using namespace cedar::cdr2;
		int k, kf, ifd;
		int updown, nstencil;

		auto & sod = const_cast<mpi::stencil_op<five_pt>&>(so);
		grid_topo & topo = sod.grid();
		MsgCtx *ctx = (MsgCtx*) sod.halo_ctx;
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		mpi::grid_func & bd = const_cast<mpi::grid_func&>(b);

		k = topo.level()+1;
		kf = topo.nlevel();
		ifd = 1;
		nstencil = 3;

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		// ibc = BMG_BCs_definite;
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		MPI_BMG2_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		                         so.len(0), so.len(1), kf, ifd, nstencil, BMG_RELAX_SYM,
		                         updown, topo.is(0), topo.is(1), ctx->msg_geom.data(),
		                         ctx->msg_geom.size(), ctx->pMSG.data(), ctx->msg_buffer.data(),
		                         ctx->msg_buffer.size(), fcomm);
	}


	template<>
	void mpi_relax_rbgs_point(const kernel_params & params,
	                          const mpi::stencil_op<nine_pt> & so,
	                          mpi::grid_func & x,
	                          const mpi::grid_func & b,
	                          const relax_stencil & sor,
	                          cycle::Dir cycle_dir)
	{
		using namespace cedar::cdr2;
		int k, kf, ifd;
		int updown, nstencil;

		auto & sod = const_cast<mpi::stencil_op<nine_pt>&>(so);
		grid_topo & topo = sod.grid();
		MsgCtx *ctx = (MsgCtx*) sod.halo_ctx;
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		mpi::grid_func & bd = const_cast<mpi::grid_func&>(b);

		k = topo.level()+1;
		kf = topo.nlevel();

		ifd = 0;
		nstencil = 5;

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		// ibc = BMG_BCs_definite;
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		MPI_BMG2_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		                         so.len(0), so.len(1), kf, ifd, nstencil, BMG_RELAX_SYM,
		                         updown, topo.is(0), topo.is(1), ctx->msg_geom.data(),
		                         ctx->msg_geom.size(), ctx->pMSG.data(), ctx->msg_buffer.data(),
		                         ctx->msg_buffer.size(), fcomm);
	}


	template<>
	void mpi_relax_lines_x(const kernel_params & params,
	                       const mpi::stencil_op<five_pt> & so,
	                       mpi::grid_func & x,
	                       const mpi::grid_func & b,
	                       const relax_stencil & sor,
	                       mpi::grid_func & res,
	                       cycle::Dir cycle_dir)
	{
		using namespace cedar::cdr2;
		int k, kf;
		int updown, nstencil;

		auto & sod = const_cast<mpi::stencil_op<five_pt>&>(so);
		grid_topo & topo = sod.grid();
		MsgCtx *ctx = (MsgCtx*) sod.halo_ctx;
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		mpi::grid_func & bd = const_cast<mpi::grid_func&>(b);

		k = topo.level()+1;
		kf = topo.nlevel();

		nstencil = 3;

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
		                              xlinecomm, ylinecomm);
	}


	template<>
	void mpi_relax_lines_x(const kernel_params & params,
	                       const mpi::stencil_op<nine_pt> & so,
	                       mpi::grid_func & x,
	                       const mpi::grid_func & b,
	                       const relax_stencil & sor,
	                       mpi::grid_func & res,
	                       cycle::Dir cycle_dir)
	{
		using namespace cedar::cdr2;
		int k, kf;
		int updown, nstencil;

		auto & sod = const_cast<mpi::stencil_op<nine_pt>&>(so);
		grid_topo & topo = sod.grid();
		MsgCtx *ctx = (MsgCtx*) sod.halo_ctx;
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		mpi::grid_func & bd = const_cast<mpi::grid_func&>(b);

		k = topo.level()+1;
		kf = topo.nlevel();

		nstencil = 5;

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
		                              xlinecomm, ylinecomm);
	}


	template<>
	void mpi_relax_lines_y(const kernel_params & params,
	                       const mpi::stencil_op<five_pt> & so,
	                       mpi::grid_func & x,
	                       const mpi::grid_func & b,
	                       const relax_stencil & sor,
	                       mpi::grid_func & res,
	                       cycle::Dir cycle_dir)
	{
		using namespace cedar::cdr2;
		int k, kf;
		int updown, nstencil;

		auto & sod = const_cast<mpi::stencil_op<five_pt>&>(so);
		grid_topo & topo = sod.grid();
		MsgCtx *ctx = (MsgCtx*) sod.halo_ctx;
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		mpi::grid_func & bd = const_cast<mpi::grid_func&>(b);

		k = topo.level()+1;
		kf = topo.nlevel();

		nstencil = 3;

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
		                              xlinecomm, ylinecomm);
	}


	template<>
	void mpi_relax_lines_y(const kernel_params & params,
	                       const mpi::stencil_op<nine_pt> & so,
	                       mpi::grid_func & x,
	                       const mpi::grid_func & b,
	                       const relax_stencil & sor,
	                       mpi::grid_func & res,
	                       cycle::Dir cycle_dir)
	{
		using namespace cedar::cdr2;
		int k, kf;
		int updown, nstencil;

		auto & sod = const_cast<mpi::stencil_op<nine_pt>&>(so);
		grid_topo & topo = sod.grid();
		MsgCtx *ctx = (MsgCtx*) sod.halo_ctx;
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		mpi::grid_func & bd = const_cast<mpi::grid_func&>(b);

		k = topo.level()+1;
		kf = topo.nlevel();

		nstencil = 5;

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
		                              xlinecomm, ylinecomm);
	}
}

}}}
