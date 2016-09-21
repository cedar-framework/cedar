#include "boxmg/2d/ftn/BMG_parameters_c.h"
#include "boxmg/2d/ftn/mpi/BMG_workspace_c.h"
#include <boxmg/3d/mpi/halo.h>


#include <boxmg/3d/cg/solve_cg.h>


extern "C" {
	using namespace boxmg;
	void BMG3_SymStd_SOLVE_cg(real_t *q, real_t *qf,
	                          len_t ii, len_t jj, len_t kk,
	                          real_t *abd, real_t *bbd, len_t nabd1, len_t nabd2,
	                          int ibc);
	void BMG3_SymStd_SOLVE_cg_LU(real_t *q, real_t *qf,
	                             len_t ii, len_t jj, len_t kk,
	                             real_t *abd, real_t *bbd,
	                             len_t nabd1, len_t nabd2, int nogm,
	                             int nproci, int nprocj, int nprock, int nproc, int rank,
	                             int *proc_grid, int *proc_coord, len_t *locarrsize,
	                             real_t *ws, len_t nmsgr, int mpicomm);
}

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void fortran_solve_cg(grid_func & x,
	                      const grid_func & b,
	                      const grid_func & ABD,
	                      real_t *bbd)
	{
		using namespace boxmg::bmg3;
		int ibc;

		auto & bd = const_cast<grid_func&>(b);
		auto & abd_data = const_cast<grid_func&>(ABD);

		ibc = BMG_BCs_definite;

		BMG3_SymStd_SOLVE_cg(x.data(), bd.data(), x.len(0), x.len(1), x.len(2),
		                     abd_data.data(), &bbd[0], ABD.len(0), ABD.len(1), ibc);
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
		BMG3_SymStd_SOLVE_cg_LU(x_par.data(), b_par.data(),
		                        x_par.len(0), x_par.len(1), x_par.len(2),
		                        abd_data.data(), bbd, abd_data.len(0), abd_data.len(1),
		                        topo.nlevel(),
		                        topo.nproc(0), topo.nproc(1), topo.nproc(2), topo.nproc(),
		                        rank, ctx->proc_grid.data(), ctx->proc_coord.data(),
		                        &ctx->msg_geom.data()[local_arr_ptr],
		                        ctx->msg_buffer.data(), ctx->msg_buffer.size(), fcomm);
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

