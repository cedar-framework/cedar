#include <mpi.h>
#include "halo.h"
#include "core/mpi/grid_func.h"
#include "fortran/BMG_parameters_c.h"
#include "fortran/mpi/BMG_workspace_c.h"

#include "solve_cg.h"


extern "C" {
	using namespace boxmg;
	void BMG2_SymStd_SOLVE_cg(real_t*, real_t*, len_t, len_t, real_t*, real_t*, len_t, len_t, int);
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
	void fortran_solve_cg(core::GridFunc & x,
	                      const core::GridFunc & b,
	                      const core::GridFunc & ABD)
	{
		using namespace boxmg::bmg2d::core;

		int ibc;

		core::GridFunc & bd = const_cast<core::GridFunc&>(b);
		core::GridFunc & abd_data = const_cast<core::GridFunc&>(ABD);


		ibc = BMG_BCs_definite;
		real_t *bbd = new real_t[ABD.len(1)];

		BMG2_SymStd_SOLVE_cg(x.data(), bd.data(), x.len(0), x.len(1),
		                    abd_data.data(), &bbd[0], ABD.len(0), ABD.len(1), ibc);

		delete[] bbd;
	}


	void solve_cg_boxmg(const solver::BoxMG &cg_solver,
	                    core::GridFunc &x,
	                    const core::GridFunc &b)
	{
		int rank;

		core::mpi::GridFunc & b_par = const_cast<core::mpi::GridFunc&>(dynamic_cast<const core::mpi::GridFunc&>(b));
		core::mpi::GridFunc & x_par = dynamic_cast<core::mpi::GridFunc&>(x);
		auto &coarse_solver = const_cast<solver::BoxMG&>(cg_solver);

		core::mpi::GridTopo & topo = b_par.grid();
		MsgCtx *ctx = (MsgCtx*) b_par.halo_ctx;


		MPI_Comm_rank(topo.comm, &rank);
		rank++; // 1 based indexing

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		len_t local_arr_ptr = ctx->pMSG(ipL_MSG_LocalArraySize,0) - 1;  // 1 vs 0 based indexing

		core::GridFunc bser(topo.nglobal(0) - 2, topo.nglobal(1) - 2);
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

		// log::set_header_msg(" (serial)");
		auto tmp = log::lvl();
		log::lvl() = 0;
		auto x_ser = coarse_solver.solve(bser);
		log::lvl() = tmp;
		//log::set_header_msg("");
		BMG2_SymStd_SOLVE_cg_unpack(x_par.data(), x_ser.data(),
		                            b.len(0), b.len(1),
		                            topo.nglobal(0), topo.nglobal(1),
		                            topo.is(0), topo.is(1),
		                            topo.nproc(0), topo.nproc(1),
		                            topo.nproc(), rank,
		                            ctx->proc_coord.data(), fcomm);

	}
}

}}}
