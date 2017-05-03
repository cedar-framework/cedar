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

class redist_solver
{
public:
	using msg_ctx = cedar::cdr3::kernel::impls::MsgCtx;
	redist_solver(const stencil_op & so, std::shared_ptr<config::reader> conf, std::array<int, 3> nblock);
	void solve(const grid_func & b, grid_func & x);

protected:
	bool ser_cg;
	std::unique_ptr<solver> slv;
	std::unique_ptr<cdr3::solver> slv_ser;
	redist_comms rcomms;
	int block_id; /** id within a block */
	int block_num; /** which block */
	std::array<int,3> nblock;
	bool active;
	int recv_id;
	std::vector<int> send_ids;
	array<len_t, len_t, 1> nbx; /** # d.o.f. for each processor in my block */
	array<len_t, len_t, 1> nby; /** # d.o.f. for each processor in my block */
	array<len_t, len_t, 1> nbz; /** # d.o.f. for each processor in my block */
	MPI_Fint msg_comm;
	grid_func b_redist;
	grid_func x_redist;
	cdr3::grid_func b_redist_ser;
	cdr3::grid_func x_redist_ser;
	std::shared_ptr<grid_topo> redist_topo(const grid_topo & fine_topo, msg_ctx & ctx);
	template <typename stencil_operator>
	stencil_operator redist_operator(const stencil_op & so, topo_ptr topo)
	{
		stencil_operator v;
		log::error << "Unspoorted type" << std::endl;
		return v;
	}
	void gather_rhs(const grid_func & b);
	void scatter_sol(grid_func & x);
	template <typename target_operator> void gather_operator(const stencil_op & src, target_operator & dest)
	{
		using buf_arr = array<len_t,real_t,1>;

		auto & sten = src.stencil();
		auto & rsten = dest.stencil();
		// save general case for later
		assert(sten.five_pt() == false);

		// Pack the operator
		buf_arr sbuf(14*sten.len(0)*sten.len(1)*sten.len(2));
		int idx = 0;
		for (auto k : sten.grange(2)) {
			for (auto j : sten.grange(1)) {
				for (auto i : sten.grange(0)) {
					sbuf(idx)   = sten(i,j,k,dir::P);
					sbuf(idx+1) = sten(i,j,k,dir::PW);
					sbuf(idx+2) = sten(i,j,k,dir::PNW);
					sbuf(idx+3) = sten(i,j,k,dir::PS);
					sbuf(idx+4) = sten(i,j,k,dir::PSW);
					sbuf(idx+5) = sten(i,j,k,dir::BNE);
					sbuf(idx+6) = sten(i,j,k,dir::BN);
					sbuf(idx+7) = sten(i,j,k,dir::BNW);
					sbuf(idx+8) = sten(i,j,k,dir::BE);
					sbuf(idx+9) = sten(i,j,k,dir::B);
					sbuf(idx+10) = sten(i,j,k,dir::BW);
					sbuf(idx+11) = sten(i,j,k,dir::BSE);
					sbuf(idx+12) = sten(i,j,k,dir::BS);
					sbuf(idx+13) = sten(i,j,k,dir::BSW);
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
								rsten(igs+ii-1,jgs+jj-1,kgs+kk-1,dir::P) = rbuf(idx);
								rsten(igs+ii-1,jgs+jj-1,kgs+kk-1,dir::PW) = rbuf(idx+1);
								rsten(igs+ii-1,jgs+jj-1,kgs+kk-1,dir::PNW) = rbuf(idx+2);
								rsten(igs+ii-1,jgs+jj-1,kgs+kk-1,dir::PS) = rbuf(idx+3);
								rsten(igs+ii-1,jgs+jj-1,kgs+kk-1,dir::PSW) = rbuf(idx+4);
								rsten(igs+ii-1,jgs+jj-1,kgs+kk-1,dir::BNE) = rbuf(idx+5);
								rsten(igs+ii-1,jgs+jj-1,kgs+kk-1,dir::BN) = rbuf(idx+6);
								rsten(igs+ii-1,jgs+jj-1,kgs+kk-1,dir::BNW) = rbuf(idx+7);
								rsten(igs+ii-1,jgs+jj-1,kgs+kk-1,dir::BE) = rbuf(idx+8);
								rsten(igs+ii-1,jgs+jj-1,kgs+kk-1,dir::B) = rbuf(idx+9);
								rsten(igs+ii-1,jgs+jj-1,kgs+kk-1,dir::BW) = rbuf(idx+10);
								rsten(igs+ii-1,jgs+jj-1,kgs+kk-1,dir::BSE) = rbuf(idx+11);
								rsten(igs+ii-1,jgs+jj-1,kgs+kk-1,dir::BS) = rbuf(idx+12);
								rsten(igs+ii-1,jgs+jj-1,kgs+kk-1,dir::BSW) = rbuf(idx+13);
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
