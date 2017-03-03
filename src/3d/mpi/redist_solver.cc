#include <cassert>
#include "boxmg/2d/ftn/mpi/BMG_workspace_c.h"

#include <boxmg/array.h>
#include <boxmg/mpi/block_partition.h>
#include <boxmg/3d/mpi/redist_solver.h>

extern "C" {
	void MSG_pause(MPI_Fint *msg_comm);
	void MSG_play(MPI_Fint msg_comm);
}

using namespace boxmg;
using namespace boxmg::bmg3::mpi;

redist_solver::redist_solver(const stencil_op & so,
                             std::shared_ptr<config::reader> conf,
                             std::array<int, 3> nblock) : nblock(nblock), active(true),
	           recv_id(-1)
{
	auto & topo = so.grid();
	msg_ctx * ctx = (msg_ctx*) so.halo_ctx;
	auto ctopo = redist_topo(topo, *ctx);
	auto rop = redist_operator(so, ctopo);
	b_redist = grid_func(ctopo);
	x_redist = grid_func(ctopo);

	if (active) {
		MPI_Fint parent_comm;
		MSG_pause(&parent_comm);
		log::push_level("redist", *conf);
		slv = std::make_unique<solver>(std::move(rop), conf);
		b_redist.halo_ctx = slv->level(-1).A.halo_ctx;
		log::pop_level();
		MSG_pause(&msg_comm);
		MSG_play(parent_comm);
	}
}


void redist_solver::solve(const grid_func & b, grid_func & x)
{
	using buf_arr = array<len_t,real_t,1>;

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
	               displs.data(), MPI_DOUBLE, collcomm);

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

	if (active) {
		MPI_Fint parent_comm;
		MSG_pause(&parent_comm);
		MSG_play(msg_comm);
		log::push_level("redist", slv->get_config());
		x_redist.set(0.0);
		slv->vcycle(x_redist, b_redist);
		log::pop_level();
		MSG_play(parent_comm);

		// copy local part from redistributed solution
		int ci = (block_id % (nbx.len(0)*nby.len(0))) % nbx.len(0);
		int cj = (block_id % (nbx.len(0)*nby.len(0))) / nbx.len(0);
		int ck = block_id / (nbx.len(0)*nby.len(0));

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

			array<len_t,real_t,3> sbuf(nbx(ci)+2, nby(cj)+2, nbz(ck)+2);
			for (auto kk : range(sbuf.len(2))) {
				for (auto jj : range(sbuf.len(1))) {
					for (auto ii : range(sbuf.len(0))) {
						sbuf(ii,jj,kk) = x_redist(igs+ii,jgs+jj,kgs+kk);
					}
				}
			}

			MPI_Send(sbuf.data(), sbuf.len(0)*sbuf.len(1)*sbuf.len(2), MPI_DOUBLE, send_id, 0, collcomm);
		}
	} else if (recv_id > -1) {
		MPI_Recv(x.data(), x.len(0)*x.len(1)*x.len(2), MPI_DOUBLE, recv_id, 0, collcomm, MPI_STATUS_IGNORE);
	}
}


stencil_op redist_solver::redist_operator(const stencil_op & so, topo_ptr topo)
{
	using buf_arr = array<len_t,real_t,1>;

	auto rop = stencil_op(topo);

	auto & sten = so.stencil();
	auto & rsten = rop.stencil();
	// save general case for later
	assert(sten.five_pt() == false);

	// Pack the operator
	buf_arr sbuf(14*sten.shape(0)*sten.shape(1)*sten.shape(2));
	int idx = 0;
	for (auto k : sten.range(2)) {
		for (auto j : sten.range(1)) {
			for (auto i : sten.range(0)) {
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
				rcounts[idx] = nbx(i)*nby(j)*nbz(k)*14;
				rbuf_len += rcounts[idx];
			}
		}
	}

	buf_arr rbuf(rbuf_len);
	MPI_Allgatherv(sbuf.data(), sbuf.len(0), MPI_DOUBLE, rbuf.data(), rcounts.data(),
	               displs.data(), MPI_DOUBLE, collcomm);

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
							rsten(igs+ii,jgs+jj,kgs+kk,dir::P) = rbuf(idx);
							rsten(igs+ii,jgs+jj,kgs+kk,dir::PW) = rbuf(idx+1);
							rsten(igs+ii,jgs+jj,kgs+kk,dir::PNW) = rbuf(idx+2);
							rsten(igs+ii,jgs+jj,kgs+kk,dir::PS) = rbuf(idx+3);
							rsten(igs+ii,jgs+jj,kgs+kk,dir::PSW) = rbuf(idx+4);
							rsten(igs+ii,jgs+jj,kgs+kk,dir::BNE) = rbuf(idx+5);
							rsten(igs+ii,jgs+jj,kgs+kk,dir::BN) = rbuf(idx+6);
							rsten(igs+ii,jgs+jj,kgs+kk,dir::BNW) = rbuf(idx+7);
							rsten(igs+ii,jgs+jj,kgs+kk,dir::BE) = rbuf(idx+8);
							rsten(igs+ii,jgs+jj,kgs+kk,dir::B) = rbuf(idx+9);
							rsten(igs+ii,jgs+jj,kgs+kk,dir::BW) = rbuf(idx+10);
							rsten(igs+ii,jgs+jj,kgs+kk,dir::BSE) = rbuf(idx+11);
							rsten(igs+ii,jgs+jj,kgs+kk,dir::BS) = rbuf(idx+12);
							rsten(igs+ii,jgs+jj,kgs+kk,dir::BSW) = rbuf(idx+13);
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

	return rop;
}


std::shared_ptr<grid_topo> redist_solver::redist_topo(const grid_topo & fine_topo, msg_ctx & ctx)
{
	using len_arr = array<int, len_t, 1>;
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

	nbx = array<len_t, len_t, 1>(high(0) - low(0) + 1);
	nby = array<len_t, len_t, 1>(high(1) - low(1) + 1);
	nbz = array<len_t, len_t, 1>(high(2) - low(2) + 1);

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
	MPI_Comm_split(fine_topo.comm, color, key, &this->collcomm);
	MPI_Comm_split(fine_topo.comm, key, color, &grid->comm);

	block_num = color;
	block_id = key;

	int nactive = (fine_topo.nproc(0) / nblock[0]) *
		(fine_topo.nproc(1) / nblock[1]) *
		(fine_topo.nproc(2) / nblock[2]);

	if (block_id > (nactive-1)) {
		active = false;
		recv_id = block_id % nactive;
	} else {
		unsigned int send_id = block_id + nactive;
		while (send_id < (nbx.len(0)*nby.len(0)*nbz.len(0))) {
			send_ids.push_back(send_id);
			send_id += nactive;
		}
	}

	return grid;
}
