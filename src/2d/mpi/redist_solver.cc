#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/mpi/redist_solver.h>

namespace cedar { namespace cdr2 { namespace mpi {

template<>
	std::unique_ptr<cdr2::stencil_op<nine_pt>> create_operator(topo_ptr topo)
{
	return std::make_unique<cdr2::stencil_op<nine_pt>>(topo->nlocal(0)-2, topo->nlocal(1)-2);
}

template<>
	std::unique_ptr<mpi::stencil_op<nine_pt>> create_operator(topo_ptr topo)
{
	return std::make_unique<mpi::stencil_op<nine_pt>>(topo);
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
template void redist_solver<cholesky_solver>::solve(grid_func & x, const grid_func & b);
template void redist_solver<multilevel_wrapper<mpi::solver<nine_pt>>>::solve(grid_func & x, const grid_func & b);
template void redist_solver<multilevel_wrapper<cdr2::solver<nine_pt>>>::solve(grid_func & x, const grid_func & b);


template<class inner_solver>
redist_solver<inner_solver>::redist_solver(const stencil_op<nine_pt> & so,
                                           service_manager<stypes> *services,
                                           std::shared_ptr<config> conf,
                                           std::array<int, 2> nblock) :
	redundant(false), nblock(nblock), active(true), recv_id(-1), services(services)
{
	// Split communicator into collective processor blocks
	auto & topo = so.grid();
	auto & halo_service = services->template get<halo_exchange>();
	auto & mpool = services->template get<mempool>();
	auto ctopo = redist_topo(topo, halo_service.leveldims(0), halo_service.leveldims(1));
	redist_operator(so, ctopo);

	if (inner_solver::is_serial) {
		b_redist = typename inner_solver::grid_func(ctopo->nlocal(0)-2, ctopo->nlocal(1)-2);
		x_redist = typename inner_solver::grid_func(ctopo->nlocal(0)-2, ctopo->nlocal(1)-2);
	} else {
		std::size_t nbytes = ctopo->nlocal(0)*ctopo->nlocal(1)*sizeof(real_t);
		real_t *xaddr = (real_t*) mpool.addr(mempool::x_redist, nbytes);
		real_t *baddr = (real_t*) mpool.addr(mempool::b_redist, nbytes);
		b_redist = grid_func(baddr, ctopo);
		x_redist = grid_func(xaddr, ctopo);
	}

	if ((redundant and active) or (not redundant and block_id == 0)) {
		MPI_Fint parent_comm;
		MSG_pause(&parent_comm);
		log::push_level("redist", *conf);
		slv = std::make_unique<inner_solver>(*so_redist, conf, *services);
		log::pop_level();
		MSG_pause(&msg_comm);
		MSG_play(parent_comm);
	}
}
template redist_solver<cholesky_solver>::redist_solver(const stencil_op<nine_pt> & so, service_manager<stypes> *services, std::shared_ptr<config> conf, std::array<int, 2> nblock);
template redist_solver<multilevel_wrapper<mpi::solver<nine_pt>>>::redist_solver(const stencil_op<nine_pt> & so, service_manager<stypes> *services, std::shared_ptr<config> conf, std::array<int, 2> nblock);
template redist_solver<multilevel_wrapper<cdr2::solver<nine_pt>>>::redist_solver(const stencil_op<nine_pt> & so, service_manager<stypes> *services, std::shared_ptr<config> conf, std::array<int, 2> nblock);


template<class inner_solver>
void redist_solver<inner_solver>::redist_operator(const stencil_op<nine_pt> & so, topo_ptr topo)
{
    auto& sod = const_cast<stencil_op<nine_pt>&>(so);
    sod.ensure_cpu();

	so_redist = create_operator<typename inner_solver::stencil_op>(topo);

	auto & rop = *so_redist;

	array<real_t,1> sbuf(5*so.len(0)*so.len(1));
	int idx = 0;
	for (auto j : so.grange(1)) {
		for (auto i : so.grange(0)) {
			sbuf(idx) = sod(i,j,nine_pt::c);
			sbuf(idx+1) = sod(i,j,nine_pt::w);
			sbuf(idx+2) = sod(i,j,nine_pt::nw);
			sbuf(idx+3) = sod(i,j,nine_pt::s);
			sbuf(idx+4) = sod(i,j,nine_pt::sw);
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
template void redist_solver<cholesky_solver>::redist_operator(const stencil_op<nine_pt> & so, topo_ptr topo);
template void redist_solver<multilevel_wrapper<mpi::solver<nine_pt>>>::redist_operator(const stencil_op<nine_pt> & so, topo_ptr topo);
template void redist_solver<multilevel_wrapper<cdr2::solver<nine_pt>>>::redist_operator(const stencil_op<nine_pt> & so, topo_ptr topo);


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

	auto & mp_service = this->services->template get<message_passing>();
	int color = grid->coord(0) + grid->nproc(0)*grid->coord(1);
	int key = (fine_topo.coord(0) - lowi) + (fine_topo.coord(1) - lowj)*parti.size(grid->coord(0));
	mp_service.comm_split(fine_topo.comm, color, key, &this->rcomms.pblock_comm);

	mp_service.comm_split(fine_topo.comm, key, color, &grid->comm);

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

	mp_service.comm_split(fine_topo.comm, color, key, &this->rcomms.active_pblock_comm);

	timer_redist(rcomms);

	return grid;
}
template std::shared_ptr<grid_topo> redist_solver<cholesky_solver>::redist_topo(const grid_topo & fine_topo, aarray<int, len_t, 2> & dimx, aarray<int, len_t, 2> & dimy);
template std::shared_ptr<grid_topo> redist_solver<multilevel_wrapper<mpi::solver<nine_pt>>>::redist_topo(const grid_topo & fine_topo, aarray<int, len_t, 2> & dimx, aarray<int, len_t, 2> & dimy);
template std::shared_ptr<grid_topo> redist_solver<multilevel_wrapper<cdr2::solver<nine_pt>>>::redist_topo(const grid_topo & fine_topo, aarray<int, len_t, 2> & dimx, aarray<int, len_t, 2> & dimy);


template<class inner_solver>
void redist_solver<inner_solver>::gather_rhs(const grid_func & b)
{
    auto& bd = const_cast<grid_func&>(b);
    bd.ensure_cpu();

	array<real_t,1> sbuf(b.shape(0)*b.shape(1));
	int idx = 0;
	for (auto j : b.range(1)) {
		for (auto i : b.range(0)) {
			sbuf(idx) = bd(i,j);
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

	auto & mp_service = this->services->template get<message_passing>();
	array<real_t,1> rbuf(rbuf_len);
	if (redundant) {
		MPI_Allgatherv(sbuf.data(), sbuf.len(0), MPI_DOUBLE, rbuf.data(), rcounts.data(),
		               displs.data(), MPI_DOUBLE, rcomms.pblock_comm);
	} else {
		mp_service.gatherv(sbuf.data(), sbuf.len(0), MPI_DOUBLE, rbuf.data(), rcounts.data(),
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
template void redist_solver<cholesky_solver>::gather_rhs(const grid_func & b);
template void redist_solver<multilevel_wrapper<mpi::solver<nine_pt>>>::gather_rhs(const grid_func & b);
template void redist_solver<multilevel_wrapper<cdr2::solver<nine_pt>>>::gather_rhs(const grid_func & b);


template<class inner_solver>
void redist_solver<inner_solver>::scatter_sol(grid_func & x)
{
    x.ensure_cpu();
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

		auto & mp_service = this->services->template get<message_passing>();
		mp_service.scatterv(sbuf.data(), scounts.data(), displs.data(),
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
        x.mark_cpu_dirty();
}
template void redist_solver<cholesky_solver>::scatter_sol(grid_func & x);
template void redist_solver<multilevel_wrapper<mpi::solver<nine_pt>>>::scatter_sol(grid_func & x);
template void redist_solver<multilevel_wrapper<cdr2::solver<nine_pt>>>::scatter_sol(grid_func & x);

}}}
