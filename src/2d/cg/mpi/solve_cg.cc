#include <mpi.h>
#include "boxmg/2d/mpi/halo.h"
#include "boxmg/2d/mpi/grid_func.h"
#include "boxmg/2d/ftn/BMG_parameters_c.h"
#include "boxmg/2d/ftn/mpi/BMG_workspace_c.h"

#include "boxmg/2d/cg/solve_cg.h"


extern "C" {
	using namespace boxmg;
	void BMG2_SymStd_SOLVE_cg_LU(real_t *Q, real_t *QF, len_t II, len_t JJ,
	                             real_t *abd, real_t *bbd,
	                             len_t nabd1, len_t nabd2, int nog,
	                             int nproci, int nprocj, int nproc, int myproc,
	                             int *proc_grid, int *proc_coord, len_t *loc_arr_size,
	                             len_t *msg_geom, len_t NMSGi, int *pMSG,
	                             real_t *msg_buffer, len_t NMSGr, int MPICOMM);
	void BMG2_SymStd_SOLVE_cg_pack(real_t *qf_ser, real_t *qf,
	                               len_t ii, len_t jj, len_t NGx, len_t NGy, len_t igs, len_t jgs,
	                               int nproci, int nprocj, int nproc, int myproc,
	                               int *proc_grid, int *proc_coord, len_t *loc_arr_size,
	                               real_t *msg_buffer, len_t nmsgr, int mpicomm);
	void BMG2_SymStd_SOLVE_cg_unpack(real_t *q, real_t *q_ser, len_t ii, len_t jj, len_t ngx,
	                                 len_t ngy, len_t igs, len_t jgs,
	                                 int nproci, int nprocj, int nproc, int myproc,
	                                 int *proc_coord, int mpicomm);
}


namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void solve_cg_boxmg(const solver &cg_solver,
	                    mpi::grid_func &x_par,
	                    const mpi::grid_func &b)
	{
		int rank;

		mpi::grid_func & b_par = const_cast<mpi::grid_func&>(b);
		auto &coarse_solver = const_cast<solver&>(cg_solver);

		grid_topo & topo = b_par.grid();
		MsgCtx *ctx = (MsgCtx*) b_par.halo_ctx;


		MPI_Comm_rank(topo.comm, &rank);
		rank++; // 1 based indexing

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		len_t local_arr_ptr = ctx->pMSG(ipL_MSG_LocalArraySize,0) - 1;  // 1 vs 0 based indexing

		grid_func bser(topo.nglobal(0) - 2, topo.nglobal(1) - 2);
		BMG2_SymStd_SOLVE_cg_pack(bser.data(), b_par.data(),
		                          b.len(0), b.len(1),
		                          topo.nglobal(0), topo.nglobal(1),
		                          topo.is(0), topo.is(1),
		                          topo.nproc(0), topo.nproc(1),
		                          topo.nproc(), rank,
		                          ctx->proc_grid.data(),
		                          ctx->proc_coord.data(),
		                          &ctx->msg_geom.data()[local_arr_ptr],
		                          ctx->msg_buffer.data(),
		                          ctx->msg_buffer.size(),
		                          fcomm);

		log::push_level("serial", coarse_solver.get_config());
		auto & x_ser = coarse_solver.level(-1).x;
		x_ser.set(0.0);
		coarse_solver.vcycle(x_ser, bser);
		log::pop_level();
		BMG2_SymStd_SOLVE_cg_unpack(x_par.data(), x_ser.data(),
		                            b.len(0), b.len(1),
		                            topo.nglobal(0), topo.nglobal(1),
		                            topo.is(0), topo.is(1),
		                            topo.nproc(0), topo.nproc(1),
		                            topo.nproc(), rank,
		                            ctx->proc_coord.data(), fcomm);

	}


	void mpi_solve_cg_lu(mpi::grid_func &x_par,
	                     const mpi::grid_func &b,
	                     const mpi::grid_func & ABD,
	                     real_t *bbd)
	{
		int rank;

		mpi::grid_func & b_par = const_cast<mpi::grid_func&>(b);
		mpi::grid_func & abd_data = const_cast<mpi::grid_func&>(ABD);

		grid_topo & topo = b_par.grid();
		MsgCtx *ctx = (MsgCtx*) b_par.halo_ctx;


		MPI_Comm_rank(topo.comm, &rank);
		rank++; // 1 based indexing

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		len_t local_arr_ptr = ctx->pMSG(ipL_MSG_LocalArraySize,0) - 1;  // 1 vs 0 based indexing
		BMG2_SymStd_SOLVE_cg_LU(x_par.data(), b_par.data(), x_par.len(0), x_par.len(1),
		                        abd_data.data(), bbd, abd_data.len(0), abd_data.len(1),
		                        topo.nlevel(),
		                        topo.nproc(0), topo.nproc(1), topo.nproc(), rank,
		                        ctx->proc_grid.data(), ctx->proc_coord.data(),
		                        &ctx->msg_geom.data()[local_arr_ptr],
		                        ctx->msg_geom.data(), ctx->msg_geom.size(),
		                        ctx->pMSG.data(), ctx->msg_buffer.data(),
		                        ctx->msg_buffer.size(), fcomm);


	}


	void solve_cg_redist(const mpi::redist_solver & cg_solver,
	                     mpi::grid_func &x,
	                     const mpi::grid_func &b)
	{
		/*
		  should move work vectors outside redist_solver
		  to eliminate need for this const_cast
		*/
		auto & slv = const_cast<mpi::redist_solver&>(cg_solver);
		slv.solve(b,x);
	}
}

}}}
