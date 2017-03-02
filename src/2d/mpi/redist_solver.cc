#include <cassert>
#include "boxmg/2d/ftn/mpi/BMG_workspace_c.h"

#include <boxmg/array.h>
#include <boxmg/mpi/block_partition.h>
#include <boxmg/2d/mpi/redist_solver.h>


extern "C" {
	void MSG_pause(MPI_Fint *msg_comm);
	void MSG_play(MPI_Fint msg_comm);
}

using namespace boxmg;
using namespace boxmg::bmg2d::mpi;

redist_solver::redist_solver(const stencil_op & so,
                             std::shared_ptr<config::reader> conf,
                             std::array<int, 2> nblock) :
	redundant(false), nblock(nblock), active(true), recv_id(-1)
{
	// Split communicator into collective processor blocks
	auto & topo = so.grid();
	msg_ctx * ctx = (msg_ctx*) so.halo_ctx;
	auto ctopo = redist_topo(topo, *ctx);
	auto rop = redist_operator(so, ctopo);
	b_redist = grid_func(ctopo);
	x_redist = grid_func(ctopo);

	if ((redundant and active) or (not redundant and block_id == 0)) {
		MPI_Fint parent_comm;
		MSG_pause(&parent_comm);
		log::push_level("redist", *conf);
		slv = std::make_unique<solver>(std::move(rop), conf);
		log::pop_level();
		MSG_pause(&msg_comm);
		MSG_play(parent_comm);
	}
}


void redist_solver::solve(const grid_func & b, grid_func & x)
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
		slv->vcycle(x_redist, b_redist);
		log::pop_level();
		MSG_play(parent_comm);
		timer_up();
	}

	timer_begin("agglomerate");
	scatter_sol(x);
	timer_end("agglomerate");
}


stencil_op redist_solver::redist_operator(const stencil_op & so, topo_ptr topo)
{
	auto rop = stencil_op(topo);

	auto & sten = so.stencil();
	auto & rsten = rop.stencil();
	// save general case for later
	assert(sten.five_pt() == false);
	rsten.five_pt() = false;

	array<len_t,real_t,1> sbuf(5*sten.shape(0)*sten.shape(1));
	int idx = 0;
	for (auto j : sten.range(1)) {
		for (auto i : sten.range(0)) {
			sbuf(idx) = sten(i,j,dir::C);
			sbuf(idx+1) = sten(i,j,dir::W);
			sbuf(idx+2) = sten(i,j,4);
			sbuf(idx+3) = sten(i,j,dir::S);
			sbuf(idx+4) = sten(i,j,dir::SW);
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
			rcounts[idx] = nbx(i)*nby(j)*5;
			rbuf_len += rcounts[idx];
		}
	}
	array<len_t,real_t,1> rbuf(rbuf_len);

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
						rsten(igs+ii,jgs+jj,dir::C) = rbuf(idx);
						rsten(igs+ii,jgs+jj,dir::W) = rbuf(idx+1);
						rsten(igs+ii,jgs+jj,4) = rbuf(idx+2);
						rsten(igs+ii,jgs+jj,dir::S) = rbuf(idx+3);
						rsten(igs+ii,jgs+jj,dir::SW) = rbuf(idx+4);
						idx += 5;
					}
				}
				igs += nx;
			}
			jgs += ny;
		}
	}

	return rop;
}


std::shared_ptr<grid_topo> redist_solver::redist_topo(const grid_topo & fine_topo, msg_ctx & ctx)
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
		grid->nlocal(0) += ctx.cg_nlocal(0, ctx.proc_grid(i, lowj)) - 2; // remove ghosts
	}
	for (auto j = lowj; j <= highj; j++) {
		grid->nlocal(1) += ctx.cg_nlocal(1, ctx.proc_grid(lowi, j)) - 2; // remove ghosts
	}
	for (unsigned int i = 0; i < lowi; i++) {
		grid->is(0) += ctx.cg_nlocal(0, ctx.proc_grid(i, lowj)) - 2; // remove ghosts
	}
	grid->is(0)++; // 1 based indexing
	for (unsigned int j = 0; j < lowj; j++) {
		grid->is(1) += ctx.cg_nlocal(1, ctx.proc_grid(lowi, j)) - 2; // remove ghosts
	}
	grid->is(1)++;

	// get ready for allgatherv
	nbx = array<len_t,len_t,1>(highi-lowi+1);
	nby = array<len_t,len_t,1>(highj-lowj+1);

	for (auto j = lowj; j <= highj; j++) {
		nby(j-lowj) = ctx.cg_nlocal(1, ctx.proc_grid(0,j)) - 2;
	}

	for (auto i = lowi; i <= highi; i++) {
		nbx(i-lowi) = ctx.cg_nlocal(0, ctx.proc_grid(i,0)) - 2;
	}

	// set dimxfine, dimyfine needed for MSG setup
	grid->dimxfine.resize(grid->nproc(0));
	for (auto i : range(grid->nproc(0))) {
		grid->dimxfine[i] = 0;
		for (auto ii = parti.low(i); ii <= parti.high(i); ii++) {
			grid->dimxfine[i] += ctx.cg_nlocal(0, ctx.proc_grid(ii, 0)) - 2;
		}
	}
	grid->dimyfine.resize(grid->nproc(1));
	for (auto j : range(grid->nproc(1))) {
		grid->dimyfine[j] = 0;
		for (auto jj = partj.low(j); jj <= partj.high(j); jj++) {
			grid->dimyfine[j] += ctx.cg_nlocal(1, ctx.proc_grid(0,jj)) - 2;
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
	int ci = fine_topo.coord(0) - lowi;
	int cj = fine_topo.coord(1) - lowj;

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


void redist_solver::gather_rhs(const grid_func & b)
{
	array<len_t,real_t,1> sbuf(b.shape(0)*b.shape(1));
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
	array<len_t,real_t,1> rbuf(rbuf_len);
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



void redist_solver::scatter_sol(grid_func & x)
{
	if (not redundant) {
		len_t sbuf_len = 0;
		for (auto j : range(nby.len(0))) {
			for (auto i : range(nbx.len(0))) {
				sbuf_len += (nbx(i)+2)*(nby(j) + 2);
			}
		}

		array<len_t, real_t, 1> sbuf(sbuf_len);
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

			array<len_t,real_t,2> sbuf(nbx(ci)+2, nby(cj)+2);
			for (auto jj : range(sbuf.len(1))) {
				for (auto ii : range(sbuf.len(0))) {
					sbuf(ii,jj) = x_redist(igs+ii,jgs+jj);
				}
			}

			MPI_Send(sbuf.data(), sbuf.len(0)*sbuf.len(1), MPI_DOUBLE, send_id, 0, rcomms.pblock_comm);
		}

	} else if (redundant and (recv_id > -1)) {
		int ci = block_id % nbx.len(0);
		int cj = block_id / nbx.len(0);

		MPI_Recv(x.data(), x.len(0)*x.len(1), MPI_DOUBLE, recv_id, 0, rcomms.pblock_comm, MPI_STATUS_IGNORE);
	}
}
