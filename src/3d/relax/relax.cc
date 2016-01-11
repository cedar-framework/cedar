#include "boxmg/2d/ftn/BMG_parameters_c.h"
#include <boxmg/3d/mpi/halo.h>

#include <boxmg/3d/relax/relax.h>

extern "C" {
	using namespace boxmg;
	void BMG3_SymStd_relax_GS(int kg, real_t *so, real_t *qf, real_t *q, real_t *sor,
	                          len_t ii, len_t jj, len_t kk, int nog, int ifd, int nstncl, int nsorv,
	                          int irelax_sym, int updown, int jpn);
	void MPI_BMG3_SymStd_relax_GS(int kg, real_t *so, real_t *qf, real_t *q, real_t *sor,
	                              len_t nlx, len_t nly, len_t nlz, len_t ngx, len_t ngy, len_t ngz,
	                              int nog, int nogm, int ifd, int nstencil, int nsorv,
	                              int irelax_sym, int updown,
	                              len_t igs, len_t jgs, len_t kgs,
	                              int myproci, int myprocj, int myprock, int myproc,
	                              int *proc_grid, int nproci, int nprocj, int nprock, int nproc,
	                              len_t *iwork, len_t nmsgi, int *pmsg, real_t *msg_buffer,
	                              len_t nmsgr, int mpicomm);
}


namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void relax_rbgs_point(const stencil_op & so,
	                      grid_func &x,
	                      const grid_func &b,
	                      const relax_stencil & sor,
	                      cycle::Dir cycle_dir)
	{
		using namespace boxmg::bmg3;
		int k, nog, ifd, nstencil, nsorv, jpn, updown;

		const grid_stencil &so_sten = so.stencil();
		auto & sod = const_cast<stencil_op&>(so);
		auto & sord = const_cast<relax_stencil&>(sor);
		auto & bd = const_cast<grid_func&>(b);

		nsorv = 2;
		// give these dummy values, ifd will pass the relevant information.
		k = 1;

		if (so_sten.five_pt()) {
			ifd = 1;
			nog = 2;
			nstencil = 4;
		} else {
			ifd = 0;
			nog = 1;
			nstencil = 14;
		}

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		jpn = BMG_BCs_definite;

		BMG3_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		                     so_sten.len(0), so_sten.len(1), so_sten.len(2), nog, ifd, nstencil, nsorv,
		                     BMG_RELAX_SYM, updown, jpn);
	}


	void mpi_relax_rbgs_point(const mpi::stencil_op & so,
	                          mpi::grid_func & x,
	                          const mpi::grid_func & b,
	                          const relax_stencil & sor,
	                          cycle::Dir cycle_dir)
	{
		using namespace boxmg::bmg3;
		int k, kf, ifd;
		int updown, nstencil;
		int rank;

		mpi::stencil_op & sod = const_cast<mpi::stencil_op&>(so);
		grid_topo & topo = sod.grid();
		MsgCtx *ctx = (MsgCtx*) sod.halo_ctx;
		grid_stencil & sten = sod.stencil();
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		mpi::grid_func & bd = const_cast<mpi::grid_func&>(b);

		k = topo.level()+1;
		kf = topo.nlevel();
		if (sten.five_pt()) {
			ifd = 1;
			nstencil = 4;
		} else {
			ifd = 0;
			nstencil = 14;
		}

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		int nsorv = 2;

		// ibc = BMG_BCs_definite;
		MPI_Comm_rank(topo.comm, &rank); rank++;
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		MPI_BMG3_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		                         sten.len(0), sten.len(1), sten.len(2),
		                         topo.nglobal(0), topo.nglobal(1), topo.nglobal(2),
		                         topo.nlevel(), topo.nlevel(), ifd, nstencil, nsorv,
		                         BMG_RELAX_SYM, updown,
		                         topo.is(0), topo.is(1), topo.is(2),
		                         topo.coord(0), topo.coord(1), topo.coord(2),
		                         rank, ctx->proc_grid.data(), topo.nproc(0), topo.nproc(1), topo.nproc(2),
		                         topo.nproc(), ctx->msg_geom.data(), ctx->msg_geom.size(), ctx->pMSG.data(),
		                         ctx->msg_buffer.data(), ctx->msg_buffer.size(), fcomm);
	}
}

}}}


