#ifndef CEDAR_2D_MPI_REDIST_SOLVER_H
#define CEDAR_2D_MPI_REDIST_SOLVER_H

#include <mpi.h>
#include <array>
#include <cassert>

#include <cedar/config.h>
#include <cedar/mpi/redist_comms.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/msg_exchanger.h>
#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"
#include <cedar/array.h>
#include <cedar/mpi/block_partition.h>
#include <cedar/util/timer.h>


extern "C" {
	void MSG_pause(MPI_Fint *msg_comm);
	void MSG_play(MPI_Fint msg_comm);
}


namespace cedar { namespace cdr2 { namespace mpi {

/**
 Coarse-grid redistribution solver.

 Provides coarse-grid redistribution through a generic solver
 interface.  This works by wrapping a solver object distributed on
 fewer cores than the callers current data distribution.

 @note This class does not currently support five_pt operators.
 */
template<class inner_solver>
class redist_solver
{
public:
	/**
	   Performs redistribution of the grid topology and stencil
	   operator given the destination 2D distribution.  This also
	   runs the setup phase of the redistributed solver.

	   @param[in] so The stencil operator to redistribute
	   @param[in] conf The config object to use for this solver
	   @param[in] nblock The destination 2D distribution
	*/
	redist_solver(const stencil_op<nine_pt> & so,
	              halo_exchanger_base *halof,
	              std::shared_ptr<config> conf,
	              std::array<int, 2> nblock);
	/**
	   Runs the redistributed solve phase
	   @param[in] b rhs
	   @param[out] x solution
	*/
	void solve(grid_func & x, const grid_func & b);

protected:
	bool redundant; /** Flag for whether redistribution is performed redundantly */
	std::unique_ptr<inner_solver> slv; /** Redistributed solver object. */
	redist_comms rcomms; /** Useful communicators used in redistribution */
	int block_id; /** id within a block */
	int block_num; /** which block */
	std::array<int,2> nblock; /** number of processors in each dimension for redistributed solve */
	bool active; /** Flag for whether the current processor is active in the redistribution */
	int recv_id; /** Where the current processor will be receiving data in the fixup phase */
	std::vector<int> send_ids; /** Where the current processor will send data in the fixup phase */
	array<len_t, 1> nbx; /** number of d.o.f. for each processor in my block */
	array<len_t, 1> nby; /** number of d.o.f. for each processor in my block */
	MPI_Fint msg_comm;
	typename inner_solver::grid_func b_redist; /** The redistributed rhs */
	typename inner_solver::grid_func x_redist; /** The redistributed solution */
	std::unique_ptr<typename inner_solver::stencil_op> so_redist; /** The redistributed operator */

	/** Redistributes the processor grid.

	    @param[in] fine_topo The source topology to redistribute
	    @param[in] dimx number of grid points on each processor. shape(proc_coord_x, level)
	    @param[in] dimy number of grid points on each processor. shape(proc_coord_y, level)
	    @returns The redistrubted processor grid topology
	*/
	std::shared_ptr<grid_topo> redist_topo(const grid_topo & fine_topo,
	                                       aarray<int, len_t, 2> & dimx,
	                                       aarray<int, len_t, 2> & dimy);
	/** Redistributes the operator.

	    @param[in] so Source operator to redistribute
	    @param[in] topo Redistributed processor grid topology
	*/
	void redist_operator(const stencil_op<nine_pt> & so, topo_ptr topo);
	/** Gathers the rhs

	    @param[in] b The rhs to redistribute
	*/
	void gather_rhs(const grid_func & b);
	/** Scatters the redistributed solution

	    @param[out] x Grid function to receieve the redistributed solution
	*/
	void scatter_sol(grid_func & x);
};


template<class inner_solver>
redist_solver<inner_solver>::redist_solver(const stencil_op<nine_pt> & so,
                                           halo_exchanger_base *halof,
                                           std::shared_ptr<config> conf,
                                           std::array<int, 2> nblock) :
	redundant(false), nblock(nblock), active(true), recv_id(-1)
{
	// Split communicator into collective processor blocks
	auto & topo = so.grid();
	auto ctopo = redist_topo(topo, halof->leveldims(0), halof->leveldims(1));
	redist_operator(so, ctopo);

	if (inner_solver::is_serial) {
		b_redist = typename inner_solver::grid_func(ctopo->nlocal(0)-2, ctopo->nlocal(1)-2);
		x_redist = typename inner_solver::grid_func(ctopo->nlocal(0)-2, ctopo->nlocal(1)-2);
	} else {
		b_redist = grid_func(ctopo);
		x_redist = grid_func(ctopo);
	}

	if ((redundant and active) or (not redundant and block_id == 0)) {
		MPI_Fint parent_comm;
		MSG_pause(&parent_comm);
		log::push_level("redist", *conf);
		slv = std::make_unique<inner_solver>(*so_redist, conf);
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

	if ((redundant and active) or (not redundant and block_id == 0)) {
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

template<class inner_solver>
void redist_solver<inner_solver>::redist_operator(const stencil_op<nine_pt> & so, topo_ptr topo)
{
	so_redist = create_operator<typename inner_solver::stencil_op>(topo);

	auto & rop = *so_redist;

	array<real_t,1> sbuf(5*so.len(0)*so.len(1));
	int idx = 0;
	for (auto j : so.grange(1)) {
		for (auto i : so.grange(0)) {
			sbuf(idx) = so(i,j,nine_pt::c);
			sbuf(idx+1) = so(i,j,nine_pt::w);
			sbuf(idx+2) = so(i,j,nine_pt::nw);
			sbuf(idx+3) = so(i,j,nine_pt::s);
			sbuf(idx+4) = so(i,j,nine_pt::sw);
			idx += 5;
		}
	}

	std::vector<int> rcounts(nbx.len(0)*nby.len(0));
	std::vector<int> displs(nbx.len(0)*nby.len(0));
	len_t rbuf_len = 0;
	for (auto j : range(nby.len(0))) {
		for (auto i : range(nbx.len(0))) {
			int idx = i+j*nbx.len(0);
			displs[idx] = rbuf_len;
			rcounts[idx] = (nbx(i)+2)*(nby(j)+2)*5;
			rbuf_len += rcounts[idx];
		}
	}
	array<real_t,1> rbuf(rbuf_len);

	if (redundant) {
		MPI_Allgatherv(sbuf.data(), sbuf.len(0), MPI_DOUBLE, rbuf.data(), rcounts.data(),
		               displs.data(), MPI_DOUBLE, rcomms.pblock_comm);
	} else {
		MPI_Gatherv(sbuf.data(), sbuf.len(0), MPI_DOUBLE, rbuf.data(), rcounts.data(),
		            displs.data(), MPI_DOUBLE, 0, rcomms.pblock_comm);
	}

	if (redundant or (block_id == 0)) {
		// Loop through all my blocks
		len_t igs, jgs;
		idx = 0;
		jgs = 1;
		for (auto j : range(nby.len(0))) {
			auto ny = 0;
			igs = 1;
			for (auto i : range(nbx.len(0))) {
				auto nx = nbx(i);
				ny = nby(j);
				for (auto jj: range(ny+2)) {
					for (auto ii : range(nx+2)) {
						rop(igs+ii-1,jgs+jj-1,nine_pt::c) = rbuf(idx);
						rop(igs+ii-1,jgs+jj-1,nine_pt::w) = rbuf(idx+1);
						rop(igs+ii-1,jgs+jj-1,nine_pt::nw) = rbuf(idx+2);
						rop(igs+ii-1,jgs+jj-1,nine_pt::s) = rbuf(idx+3);
						rop(igs+ii-1,jgs+jj-1,nine_pt::sw) = rbuf(idx+4);
						idx += 5;
					}
				}
				igs += nx;
			}
			jgs += ny;
		}
	}
}


template<class inner_solver>
	std::shared_ptr<grid_topo> redist_solver<inner_solver>::redist_topo(const grid_topo & fine_topo,
	                                                                    aarray<int, len_t, 2> & dimx,
	                                                                    aarray<int, len_t, 2> & dimy)
{
	// std::cout << fine_topo.coord(0) << " " << fine_topo.coord(1) << " => ("
	//           << fine_topo.nglobal(0) << ", " << fine_topo.nglobal(1) << ") ("
	//           << fine_topo.nlocal(0) << ", " << fine_topo.nlocal(1) << ")" << std::endl;
	auto igrd = std::make_shared<std::vector<len_t>>(NBMG_pIGRD);
	auto grid = std::make_shared<grid_topo>(igrd, 0, 1);

	grid->nproc(0) = nblock[0];
	grid->nproc(1) = nblock[1];
	grid->nproc(2) = 1;

	grid->nglobal(0) = fine_topo.nglobal(0);
	grid->nglobal(1) = fine_topo.nglobal(1);

	grid->nlocal(0) = 0;
	grid->nlocal(1) = 0;

	grid->is(0) = 0;
	grid->is(1) = 0;

	// block mapping
	block_partition parti(fine_topo.nproc(0), nblock[0]);
	block_partition partj(fine_topo.nproc(1), nblock[1]);

	grid->coord(0) = parti.owner(fine_topo.coord(0));
	grid->coord(1) = partj.owner(fine_topo.coord(1));

	auto lowi = parti.low(grid->coord(0));
	auto highi = parti.high(grid->coord(0));
	auto lowj = partj.low(grid->coord(1));
	auto highj = partj.high(grid->coord(1));
	for (auto i = lowi; i <= highi; i++) {
		grid->nlocal(0) += dimx(i, 0);
	}
	for (auto j = lowj; j <= highj; j++) {
		grid->nlocal(1) += dimy(j, 0);
	}
	for (unsigned int i = 0; i < lowi; i++) {
		grid->is(0) += dimx(i, 0);
	}
	grid->is(0)++; // 1 based indexing
	for (unsigned int j = 0; j < lowj; j++) {
		grid->is(1) += dimy(j, 0);
	}
	grid->is(1)++;

	// get ready for allgatherv
	nbx = array<len_t,1>(highi-lowi+1);
	nby = array<len_t,1>(highj-lowj+1);

	for (auto j = lowj; j <= highj; j++) {
		nby(j-lowj) = dimy(j, 0);
	}

	for (auto i = lowi; i <= highi; i++) {
		nbx(i-lowi) = dimx(i, 0);
	}

	// set dimxfine, dimyfine needed for MSG setup
	grid->dimxfine.resize(grid->nproc(0));
	for (auto i : range(grid->nproc(0))) {
		grid->dimxfine[i] = 0;
		for (auto ii = parti.low(i); ii <= parti.high(i); ii++) {
			grid->dimxfine[i] += dimx(ii, 0);
		}
	}
	grid->dimyfine.resize(grid->nproc(1));
	for (auto j : range(grid->nproc(1))) {
		grid->dimyfine[j] = 0;
		for (auto jj = partj.low(j); jj <= partj.high(j); jj++) {
			grid->dimyfine[j] += dimy(jj, 0);
		}
	}

	// add ghosts
	grid->nlocal(0) += 2;
	grid->nlocal(1) += 2;

	int color = grid->coord(0) + grid->nproc(0)*grid->coord(1);
	int key = (fine_topo.coord(0) - lowi) + (fine_topo.coord(1) - lowj)*parti.size(grid->coord(0));
	MPI_Comm_split(fine_topo.comm, color, key, &this->rcomms.pblock_comm);

	MPI_Comm_split(fine_topo.comm, key, color, &grid->comm);

	rcomms.redist_comm = grid->comm;
	rcomms.parent_comm = fine_topo.comm;

	block_num = color;
	block_id = key;

	// Find out if my processor will be included in redundant solve
	// int ci = fine_topo.coord(0) - lowi;
	// int cj = fine_topo.coord(1) - lowj;

	int nactivex = (fine_topo.nproc(0) / nblock[0]);
	int nactivey = (fine_topo.nproc(1) / nblock[1]);
	int nactive = nactivex*nactivey;
	if (block_id > (nactive-1)) {
		active = false;
		recv_id = block_id % nactive;
		color = grid->nproc(0)*grid->nproc(1);  // color for inactive processors
	} else {
		unsigned int send_id = block_id + nactive;
		while (send_id < (nbx.len(0)*nby.len(0))) {
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
	array<real_t,1> sbuf(b.shape(0)*b.shape(1));
	int idx = 0;
	for (auto j : b.range(1)) {
		for (auto i : b.range(0)) {
			sbuf(idx) = b(i,j);
			idx++;
		}
	}


	std::vector<int> rcounts(nbx.len(0)*nby.len(0));
	std::vector<int> displs(nbx.len(0)*nby.len(0));
	len_t rbuf_len = 0;
	for (auto j : range(nby.len(0))) {
		for (auto i : range(nbx.len(0))) {
			int idx = i+j*nbx.len(0);
			displs[idx] = rbuf_len;
			rcounts[idx] = nbx(i)*nby(j);
			rbuf_len += rcounts[idx];
		}
	}
	array<real_t,1> rbuf(rbuf_len);
	if (redundant) {
		MPI_Allgatherv(sbuf.data(), sbuf.len(0), MPI_DOUBLE, rbuf.data(), rcounts.data(),
		               displs.data(), MPI_DOUBLE, rcomms.pblock_comm);
	} else {
		MPI_Gatherv(sbuf.data(), sbuf.len(0), MPI_DOUBLE, rbuf.data(), rcounts.data(),
		            displs.data(), MPI_DOUBLE, 0, rcomms.pblock_comm);
	}

	if (redundant or (block_id == 0)) {
		// Loop through all my blocks
		len_t igs, jgs;
		idx = 0;
		jgs = 1;
		for (auto j : range(nby.len(0))) {
			auto ny = 0;
			igs = 1;
			for (auto i : range(nbx.len(0))) {
				auto nx = nbx(i);
				ny = nby(j);
				for (auto jj: range(ny)) {
					for (auto ii : range(nx)) {
						b_redist(igs+ii,jgs+jj) = rbuf(idx);
						idx++;
					}
				}
				igs += nx;
			}
			jgs += ny;
		}
	}
}


template<class inner_solver>
void redist_solver<inner_solver>::scatter_sol(grid_func & x)
{
	if (not redundant) {
		len_t sbuf_len = 0;
		for (auto j : range(nby.len(0))) {
			for (auto i : range(nbx.len(0))) {
				sbuf_len += (nbx(i)+2)*(nby(j) + 2);
			}
		}

		array<real_t, 1> sbuf(sbuf_len);
		std::vector<int> scounts(nbx.len(0)*nby.len(0));
		std::vector<int> displs(nbx.len(0)*nby.len(0));

		{
			len_t igs, jgs;
			jgs = 1;
			len_t idx = 0;
			for (auto j : range(nby.len(0))) {
				auto ny = 0;
				igs = 1;
				for (auto i : range(nbx.len(0))) {
					auto nx = nbx(i);
					ny = nby(j);

					int block_ind = i+j*nbx.len(0);
					displs[block_ind] = idx;
					scounts[block_ind] = (nx+2)*(ny+2);

					for (auto jj : range(ny+2)) {
						for (auto ii : range(nx+2)) {
							sbuf(idx) = x_redist(igs+ii-1,jgs+jj-1);
							idx++;
						}
					}
					igs += nx;
				}
				jgs += ny;
			}
		}
		MPI_Scatterv(sbuf.data(), scounts.data(), displs.data(),
		             MPI_DOUBLE, x.data(), x.len(0)*x.len(1),
		             MPI_DOUBLE, 0, rcomms.pblock_comm);
	}

	if (redundant and active) {

		// copy local part from redistributed solution
		int ci = block_id % nbx.len(0);
		int cj = block_id / nbx.len(0);

		len_t igs, jgs;

		igs = 1;
		for (auto i = 0; i < ci; i++) {
			igs += nbx(i);
		}
		jgs = 1;
		for (auto j = 0; j < cj; j++) {
			jgs += nby(j);
		}

		igs--; jgs--; // include ghosts

		for (auto jj : range(x.len(1))) {
			for (auto ii : range(x.len(0))) {
				x(ii,jj) = x_redist(igs+ii,jgs+jj);
			}
		}

		for (auto send_id : send_ids) {
			ci = send_id % nbx.len(0);
			cj = send_id / nbx.len(0);

			igs = 1;
			for (auto i = 0; i < ci; i++) {
				igs += nbx(i);
			}
			jgs = 1;
			for (auto j = 0; j < cj; j++) {
				jgs += nby(j);
			}

			igs--; jgs--; // include ghosts

			array<real_t,2> sbuf(nbx(ci)+2, nby(cj)+2);
			for (auto jj : range(sbuf.len(1))) {
				for (auto ii : range(sbuf.len(0))) {
					sbuf(ii,jj) = x_redist(igs+ii,jgs+jj);
				}
			}

			MPI_Send(sbuf.data(), sbuf.len(0)*sbuf.len(1), MPI_DOUBLE, send_id, 0, rcomms.pblock_comm);
		}

	} else if (redundant and (recv_id > -1)) {
		// int ci = block_id % nbx.len(0);
		// int cj = block_id / nbx.len(0);

		MPI_Recv(x.data(), x.len(0)*x.len(1), MPI_DOUBLE, recv_id, 0, rcomms.pblock_comm, MPI_STATUS_IGNORE);
	}
}

}}}

#endif
