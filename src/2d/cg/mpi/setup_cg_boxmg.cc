#include <mpi.h>
#include <cedar/2d/mpi/halo.h>
#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/types.h>

#include <cedar/2d/cg/setup_cg_boxmg.h>


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
	static void update_periodic(stencil_op<nine_pt> & o, const std::array<bool, 3> & per)
	{
		if (per[0]) {
			len_t ibeg = 1;
			len_t iend = o.shape(0);

			for (auto j : o.grange(1)) {
				o(ibeg-1,j,nine_pt::c ) = o(iend,j,nine_pt::c );
				o(ibeg-1,j,nine_pt::w ) = o(iend,j,nine_pt::w );
				o(ibeg-1,j,nine_pt::s ) = o(iend,j,nine_pt::s );
				o(ibeg-1,j,nine_pt::sw) = o(iend,j,nine_pt::sw);

				o(ibeg  ,j,nine_pt::nw) = o(iend+1,j,nine_pt::nw);
				o(ibeg-1,j,nine_pt::nw) = o(iend,j,nine_pt::nw);

				o(iend+1,j,nine_pt::c ) = o(ibeg,j,nine_pt::c );
				o(iend+1,j,nine_pt::w ) = o(ibeg,j,nine_pt::w );
				o(iend+1,j,nine_pt::s ) = o(ibeg,j,nine_pt::s );
				o(iend+1,j,nine_pt::sw) = o(ibeg,j,nine_pt::sw);
			}
		}

		if (per[1]) {
			len_t jbeg = 1;
			len_t jend = o.shape(1);

			for (auto i : o.range(0)) {
				o(i,jbeg-1,nine_pt::c ) = o(i,jend,nine_pt::c );
				o(i,jbeg-1,nine_pt::w ) = o(i,jend,nine_pt::w );
				o(i,jbeg-1,nine_pt::s ) = o(i,jend,nine_pt::s );
				o(i,jbeg-1,nine_pt::sw) = o(i,jend,nine_pt::sw);
				o(i+1,jbeg-1,nine_pt::nw) = o(i+1,jend,nine_pt::nw);

				o(i,jend+1,nine_pt::c ) = o(i,jbeg,nine_pt::c );
				o(i,jend+1,nine_pt::w ) = o(i,jbeg,nine_pt::w );
				o(i,jend+1,nine_pt::s ) = o(i,jbeg,nine_pt::s );
				o(i,jend+1,nine_pt::sw) = o(i,jbeg,nine_pt::sw);
				o(i+1,jend+1,nine_pt::nw) = o(i+1,jbeg,nine_pt::nw);
			}
		}
	}


	template<>
	void setup_cg_boxmg(const kernel_params & params,
	                    halo_exchanger * halof,
	                    const mpi::stencil_op<five_pt> & so,
	                    std::shared_ptr<config::reader> conf,
	                    std::shared_ptr<solver<five_pt>> *bmg)
	{
		log::error << "boxmg cg solver is not implemented for five point stencil" << std::endl;
	}


	template<>
	void setup_cg_boxmg(const kernel_params & params,
	                    halo_exchanger * halof,
	                    const mpi::stencil_op<nine_pt> & so,
	                    std::shared_ptr<config::reader> conf,
	                    std::shared_ptr<solver<nine_pt>> *bmg)
	{
		int nstencil, nog, rank;
		auto & sod = const_cast<mpi::stencil_op<nine_pt>&>(so);
		grid_topo & topo = sod.grid();
		MsgCtx *ctx = (MsgCtx*) halof->context_ptr();

		auto so_ser = std::make_unique<cdr2::stencil_op<nine_pt>>(topo.nglobal(0)-2, topo.nglobal(1)-2);

		nstencil = 5;

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
		                           so_ser->data(), topo.nproc(0),
		                           topo.nproc(1), topo.nproc(), rank,
		                           ctx->proc_grid.data(),
		                           ctx->proc_coord.data(),
		                           &ctx->msg_geom.data()[local_arr_ptr],
		                           fcomm);

		if (params.per_mask()) update_periodic(*so_ser, params.periodic);

		*bmg = std::make_shared<solver<nine_pt>>(*so_ser, conf);
		(*bmg)->give_op(std::move(so_ser)); // transfer ownership of this pointer
		auto & clevel = (*bmg)->levels.get(0);
		clevel.x = grid_func::like(clevel.res);
	}
}

}}}
