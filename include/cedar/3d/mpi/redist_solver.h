#ifndef CEDAR_3D_MPI_REDIST_SOLVER_H
#define CEDAR_3D_MPI_REDIST_SOLVER_H

#include <mpi.h>
#include <array>

#include <cedar/config.h>
#include <cedar/mpi/redist_comms.h>
#include <cedar/mpi/block_partition.h>
#include <cedar/3d/solver.h>
#include <cedar/3d/mpi/solver.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/mpi/msg_exchanger.h>
#include <cedar/2d/ftn/mpi/BMG_workspace_c.h>


extern "C" {
	void MSG_pause(MPI_Fint *msg_comm);
	void MSG_play(MPI_Fint msg_comm);
}


namespace cedar { namespace cdr3 { namespace mpi {

/**
 Coarse-grid redistribution solver.

 Provides coarse-grid redistribution through a generic solver
 interface.  This works by wrapping a solver object distributed on
 fewer cores than the callers current data distribution.

 @note This class does not currently support xxvii_pt operators.
 */
template<class inner_solver>
class redist_solver
{
public:
	using msg_ctx = cedar::cdr3::impls::MsgCtx;
	/**
	   Performs redistribution of the grid topology and stencil
	   operator given the destination 3D distribution.  This also
	   runs the setup phase of the redistributed solver.

	   @param[in] so The stencil operator to redistribute
	   @param[in] conf The config object to use for this solver
	   @param[in] nblock The destination 3D distribution
	*/
	redist_solver(const stencil_op<xxvii_pt> & so, mpi::msg_exchanger *halof,
	              std::shared_ptr<config> conf, std::array<int, 3> nblock);
	/**
	   Runs the redistributed solve phase
	   @param[in] b rhs
	   @param[out] x solution
	*/
	void solve(grid_func & x, const grid_func & b);

protected:
	std::unique_ptr<inner_solver> slv;
	redist_comms rcomms;
	int block_id; /** id within a block */
	int block_num; /** which block */
	std::array<int,3> nblock;
	bool active;
	int recv_id;
	std::vector<int> send_ids;
	array< len_t, 1> nbx; /** number of d.o.f. for each processor in my block */
	array< len_t, 1> nby; /** number of d.o.f. for each processor in my block */
	array< len_t, 1> nbz; /** number of d.o.f. for each processor in my block */
	MPI_Fint msg_comm;
	typename inner_solver::grid_func b_redist;
	typename inner_solver::grid_func x_redist;
	std::unique_ptr<typename inner_solver::stencil_op> so_redist; /** The redistributed operator */

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
	*/
	void redist_operator(const stencil_op<xxvii_pt> & so, topo_ptr topo);
	void gather_rhs(const grid_func & b);
	void scatter_sol(grid_func & x);
};


template<class inner_solver>
	redist_solver<inner_solver>::redist_solver(const stencil_op<xxvii_pt> & so,
	                                           mpi::msg_exchanger *halof,
	                                           std::shared_ptr<config> conf,
	                                           std::array<int, 3> nblock) :
nblock(nblock), active(true), recv_id(-1)
{
	auto & topo = so.grid();
	msg_ctx * ctx = (msg_ctx*) halof->context_ptr();
	auto ctopo = redist_topo(topo, *ctx);
	redist_operator(so, ctopo);

	if (inner_solver::is_serial) {
		b_redist = typename inner_solver::grid_func(ctopo->nlocal(0)-2, ctopo->nlocal(1)-2, ctopo->nlocal(2)-2);
		x_redist = typename inner_solver::grid_func(ctopo->nlocal(0)-2, ctopo->nlocal(1)-2, ctopo->nlocal(2)-2);
	} else {
		b_redist = grid_func(ctopo);
		x_redist = grid_func(ctopo);
	}

	if (active) {
		MPI_Fint parent_comm;
		MSG_pause(&parent_comm);
		log::push_level("redist", *conf);
		timer_down();
		slv = std::make_unique<inner_solver>(*so_redist, conf);
		timer_up();
		log::pop_level();
		MSG_pause(&msg_comm);
		MSG_play(parent_comm);
	}
}

template<class inner_solver>
	void redist_solver<inner_solver>::solve(grid_func & x, const grid_func & b)
{
	timer_begin("agglomerate");
	gather_rhs(b);
	timer_end("agglomerate");


	if (active) {
		timer_down();
		MPI_Fint parent_comm;
		MSG_pause(&parent_comm);
		MSG_play(msg_comm);
		log::push_level("redist", slv->get_config());
		x_redist.set(0.0);
		slv->cycle(x_redist, b_redist);
		log::pop_level();
		MSG_play(parent_comm);
		timer_up();
	}

	timer_begin("agglomerate");
	scatter_sol(x);
	timer_end("agglomerate");
}


template <class T>
	std::unique_ptr<T> create_operator(topo_ptr topo);
template<>
	std::unique_ptr<cdr3::stencil_op<xxvii_pt>> create_operator(topo_ptr topo)
{
	return std::make_unique<cdr3::stencil_op<xxvii_pt>>(topo->nlocal(0)-2, topo->nlocal(1)-2, topo->nlocal(2)-2);
}

template<>
	std::unique_ptr<mpi::stencil_op<xxvii_pt>> create_operator(topo_ptr topo)
{
	return std::make_unique<mpi::stencil_op<xxvii_pt>>(topo);
}


template<class inner_solver>
void redist_solver<inner_solver>::redist_operator(const stencil_op<xxvii_pt> & src,
                                                  topo_ptr topo)
{
	so_redist = create_operator<typename inner_solver::stencil_op>(topo);

	auto & dest = *so_redist;

	using buf_arr = array<real_t,1>;

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


template<class inner_solver>
	std::shared_ptr<grid_topo> redist_solver<inner_solver>::redist_topo(const grid_topo & fine_topo,
	                                                                    msg_ctx & ctx)
{
	using len_arr = aarray<int, len_t, 1>;
	auto igrd = std::make_shared<std::vector<len_t>>(NBMG_pIGRD);
	auto grid = std::make_shared<grid_topo>(igrd, 0, 1);

	block_partition part[3] = {
		block_partition(fine_topo.nproc(0), nblock[0]),
		block_partition(fine_topo.nproc(1), nblock[1]),
		block_partition(fine_topo.nproc(2), nblock[2])
	};

	len_arr low(3);
	len_arr high(3);

	for (auto i : range<int>(3)) {
		grid->nproc(i) = nblock[i];
		grid->nglobal(i) = fine_topo.nglobal(i);
		grid->nlocal(i) = 0;
		grid->is(i) = 0;
		grid->coord(i) = part[i].owner(fine_topo.coord(i));
		low(i) = part[i].low(grid->coord(i));
		high(i) = part[i].high(grid->coord(i));
	}

	for (auto i = low(0); i <= high(0); i++) {
		grid->nlocal(0) += ctx.cg_nlocal(0, ctx.proc_grid(i, low(1), low(2))) - 2; // remove ghosts
	}
	for (auto j = low(1); j <= high(1); j++) {
		grid->nlocal(1) += ctx.cg_nlocal(1, ctx.proc_grid(low(0), j, low(2))) - 2; // remove ghosts
	}
	for (auto k = low(2); k <= high(2); k++) {
		grid->nlocal(2) += ctx.cg_nlocal(2, ctx.proc_grid(low(0), low(1), k)) - 2; // remove ghosts
	}

	for (len_t i = 0; i < low(0); i++) {
		grid->is(0) += ctx.cg_nlocal(0, ctx.proc_grid(i, low(1), low(2))) - 2;
	}
	grid->is(0)++;
	for (len_t j = 0; j < low(1); j++) {
		grid->is(1) += ctx.cg_nlocal(1, ctx.proc_grid(low(0), j, low(2))) - 2;
	}
	grid->is(1)++;
	for (len_t k = 0; k < low(2); k++) {
		grid->is(2) += ctx.cg_nlocal(2, ctx.proc_grid(low(0), low(1), k)) - 2;
	}
	grid->is(2)++;

	nbx = array<len_t, 1>(high(0) - low(0) + 1);
	nby = array<len_t, 1>(high(1) - low(1) + 1);
	nbz = array<len_t, 1>(high(2) - low(2) + 1);

	for (auto i = low(0); i <= high(0); i++) {
		nbx(i-low(0)) = ctx.cg_nlocal(0, ctx.proc_grid(i,0,0)) - 2;
	}

	for (auto j = low(1); j <= high(1); j++) {
		nby(j-low(1)) = ctx.cg_nlocal(1, ctx.proc_grid(0,j,0)) - 2;
	}

	for (auto k = low(2); k <= high(2); k++) {
		nbz(k-low(2)) = ctx.cg_nlocal(2, ctx.proc_grid(0,0,k)) - 2;
	}

	grid->dimxfine.resize(grid->nproc(0));
	for (auto i : range(grid->nproc(0))) {
		grid->dimxfine[i] = 0;
		for (auto ii = part[0].low(i); ii <= part[0].high(i); ii++) {
			grid->dimxfine[i] += ctx.cg_nlocal(0, ctx.proc_grid(ii,0,0)) - 2;
		}
	}
	grid->dimyfine.resize(grid->nproc(1));
	for (auto j : range(grid->nproc(1))) {
		grid->dimyfine[j] = 0;
		for (auto jj = part[1].low(j); jj <= part[1].high(j); jj++) {
			grid->dimyfine[j] += ctx.cg_nlocal(1, ctx.proc_grid(0,jj,0)) - 2;
		}
	}
	grid->dimzfine.resize(grid->nproc(2));
	for (auto k : range(grid->nproc(2))) {
		grid->dimzfine[k] = 0;
		for (auto kk = part[2].low(k); kk <= part[2].high(k); kk++) {
			grid->dimzfine[k] += ctx.cg_nlocal(2, ctx.proc_grid(0,0,kk)) - 2;
		}
	}

	// add ghosts
	grid->nlocal(0) += 2;
	grid->nlocal(1) += 2;
	grid->nlocal(2) += 2;

	int color = grid->coord(0) + grid->nproc(0)*grid->coord(1) + grid->nproc(0)*grid->nproc(1)*grid->coord(2);
	int key = (fine_topo.coord(0) - low(0)) + (fine_topo.coord(1) - low(1))*part[0].size(grid->coord(0)) +
		(fine_topo.coord(2) - low(2)) * part[0].size(grid->coord(0))*part[1].size(grid->coord(1));
	MPI_Comm_split(fine_topo.comm, color, key, &this->rcomms.pblock_comm);
	MPI_Comm_split(fine_topo.comm, key, color, &grid->comm);

	block_num = color;
	block_id = key;

	int nactive = (fine_topo.nproc(0) / nblock[0]) *
		(fine_topo.nproc(1) / nblock[1]) *
		(fine_topo.nproc(2) / nblock[2]);

	if (block_id > (nactive-1)) {
		active = false;
		recv_id = block_id % nactive;
		color = grid->nproc(0)*grid->nproc(1)*grid->nproc(2); // color for inactive processors
	} else {
		unsigned int send_id = block_id + nactive;
		while (send_id < (nbx.len(0)*nby.len(0)*nbz.len(0))) {
			send_ids.push_back(send_id);
			send_id += nactive;
		}
	}

	MPI_Comm_split(fine_topo.comm, color, key, &this->rcomms.active_pblock_comm);

	timer_redist(rcomms);

	return grid;
}


template<class inner_solver>
void redist_solver<inner_solver>::gather_rhs(const grid_func & b)
{
	using buf_arr = array<real_t,1>;

	buf_arr sbuf(b.shape(0)*b.shape(1)*b.shape(2));
	int idx = 0;

	for (auto k : b.range(2)) {
		for (auto j : b.range(1)) {
			for (auto i : b.range(0)) {
				sbuf(idx) = b(i,j,k);
				idx++;
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
				rcounts[idx] = nbx(i)*nby(j)*nbz(k);
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
				for (auto kk : range(nz)) {
					for (auto jj : range(ny)) {
						for (auto ii : range(nx)) {
							b_redist(igs+ii,jgs+jj,kgs+kk) = rbuf(idx);
							idx++;
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


template<class inner_solver>
void redist_solver<inner_solver>::scatter_sol(grid_func & x)
{
	if (active) {
		// copy local part from redistributed solution
		int ci = (block_id % (nbx.len(0)*nby.len(0))) % nbx.len(0);
		int cj = (block_id % (nbx.len(0)*nby.len(0))) / nbx.len(0);
		int ck = block_id / (nbx.len(0)*nby.len(0));

		len_t igs, jgs, kgs;

		igs = 1;
		for (auto i = 0; i < ci; i++) {
			igs += nbx(i);
		}
		jgs = 1;
		for (auto j = 0; j < cj; j++) {
			jgs += nby(j);
		}
		kgs = 1;
		for (auto k = 0; k < ck; k++) {
			kgs += nbz(k);
		}

		igs--; jgs--; kgs--;// include ghosts

		for (auto kk : range(x.len(2))) {
			for (auto jj : range(x.len(1))) {
				for (auto ii : range(x.len(0))) {
					x(ii,jj,kk) = x_redist(igs+ii,jgs+jj,kgs+kk);
				}
			}
		}


		for (auto send_id : send_ids) {
			ci = (send_id % (nbx.len(0)*nby.len(0))) % nbx.len(0);
			cj = (send_id % (nbx.len(0)*nby.len(0))) / nbx.len(0);
			ck = send_id / (nbx.len(0)*nby.len(0));

			igs = 1;
			for (auto i = 0; i < ci; i++) {
				igs += nbx(i);
			}
			jgs = 1;
			for (auto j = 0; j < cj; j++) {
				jgs += nby(j);
			}
			kgs = 1;
			for (auto k = 0; k < ck; k++) {
				kgs += nbz(k);
			}

			igs--; jgs--; kgs--; // include ghosts

			array<real_t,3> sbuf(nbx(ci)+2, nby(cj)+2, nbz(ck)+2);
			for (auto kk : range(sbuf.len(2))) {
				for (auto jj : range(sbuf.len(1))) {
					for (auto ii : range(sbuf.len(0))) {
						sbuf(ii,jj,kk) = x_redist(igs+ii,jgs+jj,kgs+kk);
					}
				}
			}

			MPI_Send(sbuf.data(), sbuf.len(0)*sbuf.len(1)*sbuf.len(2), MPI_DOUBLE, send_id, 0, rcomms.pblock_comm);
		}
	} else if (recv_id > -1) {
		MPI_Recv(x.data(), x.len(0)*x.len(1)*x.len(2), MPI_DOUBLE, recv_id, 0, rcomms.pblock_comm, MPI_STATUS_IGNORE);
	}
}


}}}

#endif
