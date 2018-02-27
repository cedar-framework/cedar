#include "cedar/2d/ftn/BMG_parameters_c.h"
#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"
#include <cedar/3d/mpi/msg_exchanger.h>

#include <cedar/3d/cg/solve_cg.h>
#include <cedar/3d/mpi/redist_solver.h>


extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SOLVE_cg_LU(real_t *q, real_t *qf,
	                             len_t ii, len_t jj, len_t kk,
	                             real_t *abd, real_t *bbd,
	                             len_t nabd1, len_t nabd2, int nogm,
	                             int nproci, int nprocj, int nprock, int nproc, int rank,
	                             int *proc_grid, int *proc_coord, len_t *locarrsize,
	                             real_t *ws, len_t nmsgr, int mpicomm);
}

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	void mpi_solve_cg_lu(const kernel_params & params,
	                     mpi::msg_exchanger *halof,
	                     mpi::grid_func &x_par,
	                     const mpi::grid_func &b,
	                     const mpi::grid_func & ABD,
	                     real_t *bbd)
	{
		int rank;

		mpi::grid_func & b_par = const_cast<mpi::grid_func&>(b);
		mpi::grid_func & abd_data = const_cast<mpi::grid_func&>(ABD);

		grid_topo & topo = b_par.grid();
		MsgCtx *ctx = (MsgCtx*)	halof->context_ptr();


		MPI_Comm_rank(topo.comm, &rank);
		rank++; // 1 based indexing

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		len_t local_arr_ptr = ctx->pMSG(ipL_MSG_LocalArraySize,0) - 1;  // 1 vs 0 based indexing
		BMG3_SymStd_SOLVE_cg_LU(x_par.data(), b_par.data(),
		                        x_par.len(0), x_par.len(1), x_par.len(2),
		                        abd_data.data(), bbd, abd_data.len(0), abd_data.len(1),
		                        topo.nlevel(),
		                        topo.nproc(0), topo.nproc(1), topo.nproc(2), topo.nproc(),
		                        rank, ctx->proc_grid.data(), ctx->proc_coord.data(),
		                        &ctx->msg_geom.data()[local_arr_ptr],
		                        ctx->msg_buffer.data(), ctx->msg_buffer.size(), fcomm);
	}


	void solve_cg_redist(const kernel_params & params,
	                     const mpi::redist_solver & cg_solver,
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

