#include <mpi.h>
#include "cedar/2d/mpi/halo.h"
#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"
#include "cedar/2d/mpi/stencil_op.h"
#include "cedar/types.h"

#include "cedar/2d/cg/setup_cg_boxmg.h"


extern "C" {
	using namespace cedar;
	void BMG2_SymStd_SETUP_cg_boxmg(len_t ii, len_t jj, len_t ngx, len_t ngy,
	                                len_t igs, len_t jgs,
	                                real_t *so, int nstncl, int nog,
	                                real_t *ws, len_t nmsgr, real_t *so_ser,
	                                int nproci, int nprocj, int nproc, int myproc,
	                                int *proc_grid, int *proc_coord, len_t *loc_array_size,
	                                MPI_Fint MPICOMM);
}


namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	void update_periodic(stencil_op & so, const std::array<bool, 3> & per)
	{
		auto & o = so.stencil();

		if (per[0]) {
			len_t ibeg = 1;
			len_t iend = o.shape(0);

			for (auto j : o.grange(1)) {
				o(ibeg-1,j,dir::C ) = o(iend,j,dir::C );
				o(ibeg-1,j,dir::W ) = o(iend,j,dir::W );
				o(ibeg-1,j,dir::S ) = o(iend,j,dir::S );
				o(ibeg-1,j,dir::SW) = o(iend,j,dir::SW);
				o(ibeg,j,4) = o(iend+1,j,4);
				o(ibeg-1,j,4) = o(iend,j,4);

				o(iend+1,j,dir::C ) = o(ibeg,j,dir::C );
				o(iend+1,j,dir::W ) = o(ibeg,j,dir::W );
				o(iend+1,j,dir::S ) = o(ibeg,j,dir::S );
				o(iend+1,j,dir::SW) = o(ibeg,j,dir::SW);
			}
		}

		if (per[1]) {
			len_t jbeg = 1;
			len_t jend = o.shape(1);

			for (auto i : o.range(0)) {
				o(i,jbeg-1,dir::C ) = o(i,jend,dir::C );
				o(i,jbeg-1,dir::W ) = o(i,jend,dir::W );
				o(i,jbeg-1,dir::S ) = o(i,jend,dir::S );
				o(i,jbeg-1,dir::SW) = o(i,jend,dir::SW);
				o(i+1,jbeg-1,4) = o(i+1,jend,4);

				o(i,jend+1,dir::C ) = o(i,jbeg,dir::C );
				o(i,jend+1,dir::W ) = o(i,jbeg,dir::W );
				o(i,jend+1,dir::S ) = o(i,jbeg,dir::S );
				o(i,jend+1,dir::SW) = o(i,jbeg,dir::SW);
				o(i+1,jend+1,4) = o(i+1,jbeg,4);
			}
		}
	}


	void setup_cg_boxmg(const kernel_params & params,
	                    const mpi::stencil_op & so,
	                    std::shared_ptr<config::reader> conf,
	                    std::shared_ptr<solver> *bmg)
	{
		int nstencil, nog, rank;
		mpi::stencil_op & sod = const_cast<mpi::stencil_op&>(so);
		grid_topo & topo = sod.grid();
		MsgCtx *ctx = (MsgCtx*) sod.halo_ctx;
		auto & sten = sod.stencil();

		stencil_op so_ser(topo.nglobal(0)-2, topo.nglobal(1)-2);
		so_ser.stencil().five_pt() = sten.five_pt();

		if (sten.five_pt()) nstencil = 3;
		else nstencil = 5;

		nog = topo.nlevel();

		MPI_Comm_rank(topo.comm, &rank);
		rank++; // 1 based indexing

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		len_t local_arr_ptr = ctx->pMSG(ipL_MSG_LocalArraySize,0) - 1;  // 1 vs 0 based indexing

		BMG2_SymStd_SETUP_cg_boxmg(topo.nlocal(0), topo.nlocal(1),
		                           topo.nglobal(0), topo.nglobal(1),
		                           topo.is(0), topo.is(1),
		                           sod.data(), nstencil, nog,
		                           ctx->msg_buffer.data(),
		                           ctx->msg_buffer.size(),
		                           so_ser.data(), topo.nproc(0),
		                           topo.nproc(1), topo.nproc(), rank,
		                           ctx->proc_grid.data(),
		                           ctx->proc_coord.data(),
		                           &ctx->msg_geom.data()[local_arr_ptr],
		                           fcomm);

		if (params.per_mask()) update_periodic(so_ser, params.periodic);

		*bmg = std::make_shared<solver>(std::move(so_ser), conf);
		(*bmg)->level(-1).x = grid_func::like((*bmg)->level(-1).res);
	}
}

}}}
