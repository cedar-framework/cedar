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

redist_solver::redist_solver(const stencil_op & so, std::array<int, 2> nblock) : nblock(nblock)
{
	// Split communicator into collective processor blocks
	auto & topo = so.grid();
	msg_ctx * ctx = (msg_ctx*) so.halo_ctx;
	auto ctopo = redist_topo(topo, *ctx);
	auto rop = redist_operator(so, ctopo);
	b_redist = grid_func(ctopo);
	x_redist = grid_func(ctopo);

//	if (block_id == 2) {
		MPI_Fint parent_comm;
		MSG_pause(&parent_comm);
		log::set_header_msg(" (redist)");
		slv = std::make_unique<solver>(std::move(rop));
		log::set_header_msg("");
		MSG_pause(&msg_comm);
		MSG_play(parent_comm);
		// if (block_num == 1) {
		// 	std::ofstream rfile;
		// 	rfile.open("after", std::ios::out | std::ios::trunc | std::ios::binary);
		// 	rfile << (*slv).level(-1).A;
		// 	rfile.close();
		// }
//	}
	// MPI_Barrier(topo.comm);
	// MPI_Abort(topo.comm,0);
}


void redist_solver::solve(const grid_func & b, grid_func & x)
{
	array<len_t,real_t,1> sbuf(b.shape(0)*b.shape(1));
	int idx = 0;
	for (auto j : b.range(1)) {
		for (auto i : b.range(0)) {
			sbuf(idx) = b(i,j);
			idx++;
		}
	}


	std::vector<int> rcounts(nlocal.len(0)*nlocal.len(1));
	std::vector<int> displs(nlocal.len(0)*nlocal.len(1));
	len_t rbuf_len = 0;
	for (auto j : range(nlocal.len(1))) {
		for (auto i : range(nlocal.len(0))) {
			int idx = i+j*nlocal.len(0);
			displs[idx] = rbuf_len;
			rcounts[idx] = nlocal(0,i,j)*nlocal(1,i,j);
			rbuf_len += rcounts[idx];
		}
	}
	array<len_t,real_t,1> rbuf(rbuf_len);
	MPI_Allgatherv(sbuf.data(), sbuf.len(0), MPI_DOUBLE, rbuf.data(), rcounts.data(),
	               displs.data(), MPI_DOUBLE, collcomm);


	// Loop through all my blocks
	len_t igs, jgs;
	idx = 0;
	jgs = 1;
	for (auto j : range(nlocal.len(2))) {
		auto ny = 0;
		igs = 1;
		for (auto i : range(nlocal.len(1))) {
			// if (i == 1 and j == 0) {
			// 	std::cout << "igs,jgs: " << igs << " " << jgs << std::endl;
			// 	std::cout << "nlx,nly: " << nlocal(0,i,j) << " " << nlocal(1,i,j) << std::endl;
			// }
			auto nx = nlocal(0,i,j);
			ny = nlocal(1,i,j);
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

	MPI_Fint parent_comm;
	MSG_pause(&parent_comm);
	MSG_play(msg_comm);
	log::set_header_msg(" (redist)");
	slv->solve(b_redist, x_redist);
	log::set_header_msg("");
	// MSG_pause(&msg_comm);
	MSG_play(parent_comm);

	int ci = block_id % nlocal.len(1);
	int cj = block_id / nlocal.len(0);

	igs = 1;
	for (auto i = 0; i < ci; i++) {
		igs += nlocal(0,i,0);
	}
	jgs = 1;
	for (auto j = 0; j < cj; j++) {
		jgs += nlocal(1,0,j);
	}

	igs--; jgs--; // include ghosts

	for (auto jj : range(x.len(1))) {
		for (auto ii : range(x.len(0))) {
			x(ii,jj) = x_redist(igs+ii,jgs+jj);
		}
	}
}


stencil_op redist_solver::redist_operator(const stencil_op & so, topo_ptr topo)
{
	auto rop = stencil_op(topo);

	auto & sten = so.stencil();
	auto & rsten = rop.stencil();
	// save general case for later
	assert(sten.five_pt() == false);

	//if (block_num == 1) {
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

	std::vector<int> rcounts(nlocal.len(0)*nlocal.len(1));
	std::vector<int> displs(nlocal.len(0)*nlocal.len(1));
	len_t rbuf_len = 0;
	for (auto j : range(nlocal.len(1))) {
		for (auto i : range(nlocal.len(0))) {
			int idx = i+j*nlocal.len(0);
			displs[idx] = rbuf_len;
			rcounts[idx] = nlocal(0,i,j)*nlocal(1,i,j)*5;
			rbuf_len += rcounts[idx];
		}
	}
	array<len_t,real_t,1> rbuf(rbuf_len);
	MPI_Allgatherv(sbuf.data(), sbuf.len(0), MPI_DOUBLE, rbuf.data(), rcounts.data(),
	               displs.data(), MPI_DOUBLE, collcomm);

	// {// DEBUG
	// 	if (block_id == 0) {
	// 		idx = rcounts[0];
	// 		for (auto j : range(nlocal(1,1,0))) {
	// 			for (auto i : range(nlocal(0,1,0))) {
	// 				for (auto k : range(5)) {
	// 					std::cout << "after: " << rbuf(idx+k) << std::endl;
	// 				}
	// 				idx += 5;
	// 			}
	// 		}
	// 	}
	// }

	// Loop through all my blocks
	len_t igs, jgs;
	idx = 0;
	jgs = 1;
	for (auto j : range(nlocal.len(2))) {
		auto ny = 0;
		igs = 1;
		for (auto i : range(nlocal.len(1))) {
			// if (i == 1 and j == 0) {
			// 	std::cout << "igs,jgs: " << igs << " " << jgs << std::endl;
			// 	std::cout << "nlx,nly: " << nlocal(0,i,j) << " " << nlocal(1,i,j) << std::endl;
			// }
			auto nx = nlocal(0,i,j);
			ny = nlocal(1,i,j);
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

	// { // debug
	// 	if (block_num == 1) {
	// 		std::ofstream cfile;
	// 		int rank;
	// 		auto &ctopo = so.grid();
	// 		MPI_Comm_rank(ctopo.comm, &rank);
	// 		cfile.open("stencil-before-" + std::to_string(rank), std::ios::out | std::ios::trunc | std::ios::binary);

	// 		cfile << so;
	// 		cfile.close();
	// 	}

	// 	// if (block_num == 0 and block_id == 0) {
	// 	// 	std::ofstream rfile;
	// 	// 	rfile.open("after", std::ios::out | std::ios::trunc | std::ios::binary);
	// 	// 	rfile << rop;
	// 	// 	rfile.close();

	// 	// 	auto &ctopo = so.grid();
	// 	// 	auto &rtopo = rop.grid();
	// 	// 	//std::cout << "after: " << ctopo.is(0) << " " << ctopo.is(1) << " " << rtopo.is(0) << " " << rtopo.is(1) << std::endl;
	// 	// 	//std::cout << ctopo.nglobal(0) << " " << ctopo.nglobal(1) << " => " << rtopo.nglobal(0) << " " << rtopo.nglobal(1) << std::endl;
	// 	// 	//std::cout << ctopo.nlocal(0) << " " << ctopo.nlocal(1) << " => " << nlocal(0,1,0) << " " << nlocal(1,1,0) << std::endl;
	// 	// }
	// }

	// {
	// 	auto & tp = rop.grid();
	// 	auto & rsten = rop.stencil();
	// 	std::cout << tp.nglobal(0) << " vs " << rsten.len(0) << std::endl;
	// }

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
	for (auto i = 0; i < lowi; i++) {
		grid->is(0) += ctx.cg_nlocal(0, ctx.proc_grid(i, lowj)) - 2; // remove ghosts
	}
	grid->is(0)++; // 1 based indexing
	for (auto j = 0; j < lowj; j++) {
		grid->is(1) += ctx.cg_nlocal(1, ctx.proc_grid(lowi, j)) - 2; // remove ghosts
	}
	grid->is(1)++;

	// get ready for allgatherv
	nlocal = array<len_t,len_t,3>(2,(highi-lowi+1), (highj-lowj+1));
	for (auto j = lowj; j <= highj; j++) {
		for (auto i = lowi; i <= highi; i++) {
			auto nx = ctx.cg_nlocal(0, ctx.proc_grid(i,j)) - 2;
			auto ny = ctx.cg_nlocal(1, ctx.proc_grid(i,j)) - 2;
			nlocal(0,i-lowi,j-lowj) = nx;
			nlocal(1,i-lowi,j-lowj) = ny;
		}
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
	MPI_Comm_split(fine_topo.comm, color, key, &this->collcomm);

	MPI_Comm_split(fine_topo.comm, key, color, &grid->comm);

	block_num = color;
	block_id = key;

	// if (block_id == 1) {
	// 	//std::cout << block_num << std::endl;
	// 	int comm_size;
	// 	int rank;
	// 	MPI_Comm_size(grid->comm, &comm_size);
	// 	MPI_Comm_rank(grid->comm, &rank);
	// 	std::cout << "[" << rank << "] " << "comm size: " << comm_size << " vs " << grid->nproc() << std::endl;
	// }

	// if (block_num == 2) {
	// 	// for (auto j = 0; j < lowj; j++) {
	// 	// 	std::cout << "[" << j << "] " << ctx.cg_nlocal(0, 
	// 	// }
	// 	// std::cout << grid->is(0) << " " << grid->is(1) << std::endl;
	// }

	// if (fine_topo.coord(0) == 0 and fine_topo.coord(1) == 0) {
	// 	std::cout << lowi << " " << highi << std::endl;
	// 	std::cout << lowj << " " << highj << std::endl;
	// 	std::cout << ctx.cg_nlocal(0, ctx.proc_grid(1,0)) - 2 << std::endl;
	// }

	// std::cout << grid->coord(0) << " " << grid->coord(1) << " => "
	//           << grid->nlocal(0)-2 << " "  << grid->nlocal(1)-2 << std::endl;

	// std::cout << fine_topo.coord(0) << " " << fine_topo.coord(1) << " => ("
	//           << color << " " << key << ")" << std::endl;

	// MPI_Barrier(fine_topo.comm);
	// MPI_Abort(fine_topo.comm, 0);

	return grid;
}
