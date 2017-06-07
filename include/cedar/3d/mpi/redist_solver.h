#ifndef CEDAR_3D_MPI_REDIST_SOLVER_H
#define CEDAR_3D_MPI_REDIST_SOLVER_H

#include <mpi.h>
#include <array>

#include <cedar/config/reader.h>
#include <cedar/mpi/redist_comms.h>
#include <cedar/3d/solver.h>
#include <cedar/3d/mpi/solver.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/mpi/halo.h>

namespace cedar { namespace cdr3 { namespace mpi {

/**
 Coarse-grid redistribution solver.

 Provides coarse-grid redistribution through a generic solver
 interface.  This works by wrapping a solver object distributed on
 fewer cores than the callers current data distribution.

 @note This class does not currently support xxvii_pt operators.
 */
class redist_solver
{
public:
	using msg_ctx = cedar::cdr3::kernel::impls::MsgCtx;
	/**
	   Performs redistribution of the grid topology and stencil
	   operator given the destination 3D distribution.  This also
	   runs the setup phase of the redistributed solver.

	   @param[in] so The stencil operator to redistribute
	   @param[in] conf The config object to use for this solver
	   @param[in] nblock The destination 3D distribution
	*/
	redist_solver(const stencil_op<xxvii_pt> & so, std::shared_ptr<config::reader> conf, std::array<int, 3> nblock);
	/**
	   Runs the redistributed solve phase
	   @param[in] b rhs
	   @param[out] x solution
	*/
	void solve(const grid_func & b, grid_func & x);

protected:
	bool ser_cg;
	std::unique_ptr<solver> slv;
	std::unique_ptr<cdr3::solver<xxvii_pt>> slv_ser;
	redist_comms rcomms;
	int block_id; /** id within a block */
	int block_num; /** which block */
	std::array<int,3> nblock;
	bool active;
	int recv_id;
	std::vector<int> send_ids;
	array<len_t, len_t, 1> nbx; /** number of d.o.f. for each processor in my block */
	array<len_t, len_t, 1> nby; /** number of d.o.f. for each processor in my block */
	array<len_t, len_t, 1> nbz; /** number of d.o.f. for each processor in my block */
	MPI_Fint msg_comm;
	grid_func b_redist;
	grid_func x_redist;
	cdr3::grid_func b_redist_ser;
	cdr3::grid_func x_redist_ser;

	/** Redistributes the processor grid.

	    @param[in] fine_topo The source topology to redistribute
	    @param[in] msg_ctx The MSG context struct
	    @returns The redistrubted processor grid topology
	*/
	std::shared_ptr<grid_topo> redist_topo(const grid_topo & fine_topo, msg_ctx & ctx);

	/** Redistributes the operator.

	    @tparam stencil_operator Type of operator to redistribute (serial or MPI)
	    @param[in] so Source operator to redistribute
	    @param[in] topo Redistributed processor grid topology
	    @returns The redistributed operator
	*/
	template <stencil_operator>
	stencil_operator redist_operator(const stencil_op<sten> & so, topo_ptr topo)
	{
		stencil_operator v;
		log::error << "Unspoorted type" << std::endl;
		return v;
	}
	void gather_rhs(const grid_func & b);
	void scatter_sol(grid_func & x);
	template <typename target_operator> void gather_operator(const stencil_op<xxvii_pt> & src,
	                                                         target_operator & dest)
	{
		using buf_arr = array<len_t,real_t,1>;

		// Pack the operator
		buf_arr sbuf(14*src.len(0)*src.len(1)*src.len(2));
		int idx = 0;
		for (auto k : src.grange(2)) {
			for (auto j : src.grange(1)) {
				for (auto i : src.grange(0)) {
					sbuf(idx)   = src(i,j,k,xxvii_pt::p  );
					sbuf(idx+1) = src(i,j,k,xxvii_pt::pw );
					sbuf(idx+2) = src(i,j,k,xxvii_pt::pnw);
					sbuf(idx+3) = src(i,j,k,xxvii_pt::ps );
					sbuf(idx+4) = src(i,j,k,xxvii_pt::psw);
					sbuf(idx+5) = src(i,j,k,xxvii_pt::bne);
					sbuf(idx+6) = src(i,j,k,xxvii_pt::bn );
					sbuf(idx+7) = src(i,j,k,xxvii_pt::bnw);
					sbuf(idx+8) = src(i,j,k,xxvii_pt::be );
					sbuf(idx+9) = src(i,j,k,xxvii_pt::b  );
					sbuf(idx+10) = src(i,j,k,xxvii_pt::bw );
					sbuf(idx+11) = src(i,j,k,xxvii_pt::bse);
					sbuf(idx+12) = src(i,j,k,xxvii_pt::bs );
					sbuf(idx+13) = src(i,j,k,xxvii_pt::bsw);
					idx += 14;
				}
			}
		}

		std::vector<int> rcounts(nbx.len(0)*nby.len(0)*nbz.len(0));
		std::vector<int> displs(nbx.len(0)*nby.len(0)*nbz.len(0));
		len_t rbuf_len = 0;
		for (auto k : range(nbz.len(0))) {
			for (auto j : range(nby.len(0))) {
				for (auto i : range(nbx.len(0))) {
					int idx = i + j*nbx.len(0) + k*nbx.len(0)*nby.len(0);
					displs[idx] = rbuf_len;
					rcounts[idx] = (nbx(i)+2)*(nby(j)+2)*(nbz(k)+2)*14;
					rbuf_len += rcounts[idx];
				}
			}
		}

		buf_arr rbuf(rbuf_len);
		MPI_Allgatherv(sbuf.data(), sbuf.len(0), MPI_DOUBLE, rbuf.data(), rcounts.data(),
		               displs.data(), MPI_DOUBLE, rcomms.pblock_comm);

		// Loop through all my blocks
		// TODO: this is unreadable, reduce the number of nestings
		len_t igs, jgs, kgs;
		idx = 0;
		kgs = 1;
		for (auto k : range(nbz.len(0))) {
			auto nz = 0;
			jgs = 1;
			for (auto j : range(nby.len(0))) {
				auto ny = 0;
				igs = 1;
				for (auto i : range(nbx.len(0))) {
					auto nx = nbx(i);
					ny = nby(j);
					nz = nbz(k);
					for (auto kk : range(nz+2)) {
						for (auto jj : range(ny+2)) {
							for (auto ii : range(nx+2)) {
								dest(igs+ii-1,jgs+jj-1,kgs+kk-1,xxvii_pt::p  ) = rbuf(idx);
								dest(igs+ii-1,jgs+jj-1,kgs+kk-1,xxvii_pt::pw ) = rbuf(idx+1);
								dest(igs+ii-1,jgs+jj-1,kgs+kk-1,xxvii_pt::pnw) = rbuf(idx+2);
								dest(igs+ii-1,jgs+jj-1,kgs+kk-1,xxvii_pt::ps ) = rbuf(idx+3);
								dest(igs+ii-1,jgs+jj-1,kgs+kk-1,xxvii_pt::psw) = rbuf(idx+4);
								dest(igs+ii-1,jgs+jj-1,kgs+kk-1,xxvii_pt::bne) = rbuf(idx+5);
								dest(igs+ii-1,jgs+jj-1,kgs+kk-1,xxvii_pt::bn ) = rbuf(idx+6);
								dest(igs+ii-1,jgs+jj-1,kgs+kk-1,xxvii_pt::bnw) = rbuf(idx+7);
								dest(igs+ii-1,jgs+jj-1,kgs+kk-1,xxvii_pt::be ) = rbuf(idx+8);
								dest(igs+ii-1,jgs+jj-1,kgs+kk-1,xxvii_pt::b  ) = rbuf(idx+9);
								dest(igs+ii-1,jgs+jj-1,kgs+kk-1,xxvii_pt::bw ) = rbuf(idx+10);
								dest(igs+ii-1,jgs+jj-1,kgs+kk-1,xxvii_pt::bse) = rbuf(idx+11);
								dest(igs+ii-1,jgs+jj-1,kgs+kk-1,xxvii_pt::bs ) = rbuf(idx+12);
								dest(igs+ii-1,jgs+jj-1,kgs+kk-1,xxvii_pt::bsw) = rbuf(idx+13);
								idx += 14;
							}
						}
					}
					igs += nx;
				}
				jgs += ny;
			}
			kgs += nz;
		}
	}
};

}}}

#endif
