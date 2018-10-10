#include <cedar/3d/mpi/plane_mempool.h>
#include <cedar/3d/mpi/plane_mpi.h>
#include <cedar/3d/mpi/plane_exchange.h>
#include <cedar/3d/mpi/relax_planes.h>


namespace cedar { namespace cdr3 { namespace mpi {

std::tuple<int, MPI_Comm> log_begin(bool log_planes, int ipl, const std::string & suff, MPI_Comm comm)
{
	auto tmp = std::make_tuple<int, MPI_Comm>(log::lvl(), log::get_comm());

	if (log_planes) {
		log::set_header_msg(" (plane-" + suff + " " + std::to_string(ipl) + ")");
		log::set_comm(comm);
	} else
		log::lvl() = 0;

	return tmp;
}


void log_end(bool log_planes, std::tuple<int, MPI_Comm> saved)
{
	if (log_planes) {
		log::set_header_msg("");
		log::set_comm(std::get<1>(saved));
	} else
		log::lvl() = std::get<0>(saved);
}


cdr2::mpi::kman_ptr master_kman(config & conf, int nplanes, bool aggregate, plane_team & team)
{
	auto kman = cdr2::mpi::build_kernel_manager(conf);
	auto & serv = kman->services();

	auto regis = [nplanes, aggregate, &team](service_manager<cdr2::mpi::stypes> & sman) {
		team.masters.push_back(&sman);
		sman.add<services::message_passing, plane_setup_mpi>("plane-setup");
		sman.add<services::mempool, plane_mempool>("plane", nplanes);
		#ifdef PLANE_AGG
		sman.add<services::message_passing, plane_mpi>("plane", nplanes, team.threads);
		sman.add<services::halo_exchange<cdr2::mpi::stypes>, plane_exchange>("plane", nplanes, team.threads);
		#endif

		sman.set<services::message_passing>("plane-setup");
		if (aggregate)
			sman.set<services::mempool>("plane");
	};

	serv.set_user_reg(regis);
	return kman;
}


cdr2::mpi::kman_ptr worker_kman(config & conf, int nplanes, bool aggregate, plane_team & team, int worker_id)
{
	auto kman = cdr2::mpi::build_kernel_manager(conf);
	auto & serv = kman->services();

	auto regis = [worker_id, nplanes, aggregate, &team](service_manager<cdr2::mpi::stypes> & sman) {
		service_manager<cdr2::mpi::stypes> & master = team.get_master();

		// add worker to team
		auto *sman_ptr = &sman;
		team.add_worker(sman_ptr);

		auto mp_service = master.get_ptr<services::message_passing>();
		auto mempool_service = master.get_ptr<services::mempool>();
		auto *mpi_keys = static_cast<plane_setup_mpi*>(mp_service.get())->get_keys();
		auto *addrs = static_cast<plane_mempool*>(mempool_service.get())->get_addrs();

		sman.add<services::message_passing, plane_setup_mpi>("plane-setup", mpi_keys);
		sman.add<services::mempool, plane_mempool>("plane", worker_id, addrs);
		#ifdef PLANE_AGG
		sman.add<services::message_passing, plane_mpi>("plane", nplanes, worker_id, team.threads);
		sman.add<services::halo_exchange<cdr2::mpi::stypes>, plane_exchange>("plane", nplanes, worker_id, team.threads);
		#endif

		sman.set<services::message_passing>("plane-setup");
		if (aggregate)
			sman.set<services::mempool>("plane");
	};

	serv.set_user_reg(regis);
	return kman;
}


template<>
void copy_rhs<relax_dir::xy>(const stencil_op<seven_pt> & so,
                             const grid_func & x,
                             const grid_func & b,
                             cdr2::mpi::grid_func & b2,
                             int ipl)
{
	for (auto j : b.range(1)) {
		for (auto i : b.range(0)) {
			b2(i,j) = b(i,j,ipl)
				+ so(i,j,ipl,seven_pt::b) * x(i,j,ipl-1)
				+ so(i,j,ipl+1,seven_pt::b) * x(i,j,ipl+1);
		}
	}
}


template<>
void copy_rhs<relax_dir::xy>(const stencil_op<xxvii_pt> & so,
                             const grid_func & x,
                             const grid_func & b,
                             cdr2::mpi::grid_func & b2,
                             int ipl)
{
	for (auto j : b.range(1)) {
		for (auto i : b.range(0)) {
			b2(i,j) = b(i,j,ipl)
				+ so(i,j,ipl,xxvii_pt::b)*x(i,j,ipl-1)
				+ so(i,j,ipl,xxvii_pt::bw)*x(i-1,j,ipl-1)
				+ so(i,j+1,ipl,xxvii_pt::bnw)*x(i-1,j+1,ipl-1)
				+ so(i,j+1,ipl,xxvii_pt::bn)*x(i,j+1,ipl-1)
				+ so(i+1,j+1,ipl,xxvii_pt::bne)*x(i+1,j+1,ipl-1)
				+ so(i+1,j,ipl,xxvii_pt::be)*x(i+1,j,ipl-1)
				+ so(i+1,j,ipl,xxvii_pt::bse)*x(i+1,j-1,ipl-1)
				+ so(i,j,ipl,xxvii_pt::bs)*x(i,j-1,ipl-1)
				+ so(i,j,ipl,xxvii_pt::bsw)*x(i-1,j-1,ipl-1)
				+ so(i,j,ipl+1,xxvii_pt::be)*x(i-1,j,ipl+1)
				+ so(i,j+1,ipl+1,xxvii_pt::bse)*x(i-1,j+1,ipl+1)
				+ so(i,j+1,ipl+1,xxvii_pt::bs)*x(i,j+1,ipl+1)
				+ so(i+1,j+1,ipl+1,xxvii_pt::bsw)*x(i+1,j+1,ipl+1)
				+ so(i+1,j,ipl+1,xxvii_pt::bw)*x(i+1,j,ipl+1)
				+ so(i,j,ipl+1,xxvii_pt::b)*x(i,j,ipl+1)
				+ so(i+1,j,ipl+1,xxvii_pt::bnw)*x(i+1,j-1,ipl+1)
				+ so(i,j,ipl+1,xxvii_pt::bn)*x(i,j-1,ipl+1)
				+ so(i,j,ipl+1,xxvii_pt::bne)*x(i-1,j-1,ipl+1);
		}
	}
}


template<>
void copy_rhs<relax_dir::xz>(const stencil_op<seven_pt> & so,
                             const grid_func & x,
                             const grid_func & b,
                             cdr2::mpi::grid_func & b2,
                             int ipl)
{
	for (auto k : b.range(2)) {
		for (auto i : b.range(0)) {
			b2(i,k) = b(i,ipl,k)
				+ so(i,ipl,k,seven_pt::ps) * x(i,ipl-1,k)
				+ so(i,ipl+1,k,seven_pt::ps) * x(i,ipl+1,k);
		}
	}
}


template<>
void copy_rhs<relax_dir::xz>(const stencil_op<xxvii_pt> & so,
                             const grid_func & x,
                             const grid_func & b,
                             cdr2::mpi::grid_func & b2,
                             int ipl)
{
	for (auto k : b.range(2)) {
		for (auto i : b.range(0)) {
			b2(i,k) = b(i,ipl,k)
				+ so(i,ipl+1,k,xxvii_pt::pnw) * x(i-1,ipl+1,k)
				+ so(i,ipl+1,k,xxvii_pt::ps) * x(i,ipl+1,k)
				+ so(i+1,ipl+1,k,xxvii_pt::psw) * x(i+1,ipl+1,k)
				+ so(i,ipl+1,k,xxvii_pt::bnw) * x(i-1,ipl+1,k-1)
				+ so(i,ipl+1,k,xxvii_pt::bn) * x(i,ipl+1,k-1)
				+ so(i+1,ipl+1,k,xxvii_pt::bne) * x(i+1,ipl+1,k-1)
				+ so(i,ipl+1,k+1,xxvii_pt::bse) * x(i-1,ipl+1,k+1)
				+ so(i,ipl+1,k+1,xxvii_pt::bs) * x(i,ipl+1,k+1)
				+ so(i+1,ipl+1,k+1,xxvii_pt::bsw) * x(i+1,ipl+1,k+1)
				+ so(i,ipl,k,xxvii_pt::psw) * x(i-1,ipl-1,k)
				+ so(i,ipl,k,xxvii_pt::ps) * x(i,ipl-1,k)
				+ so(i+1,ipl,k,xxvii_pt::pnw) * x(i+1,ipl-1,k)
				+ so(i,ipl,k,xxvii_pt::bsw) * x(i-1,ipl-1,k-1)
				+ so(i,ipl,k,xxvii_pt::bs) * x(i,ipl-1,k-1)
				+ so(i+1,ipl,k,xxvii_pt::bse) * x(i+1,ipl-1,k-1)
				+ so(i,ipl,k+1,xxvii_pt::bne) * x(i-1,ipl-1,k+1)
				+ so(i,ipl,k+1,xxvii_pt::bn) * x(i,ipl-1,k+1)
				+ so(i+1,ipl,k+1,xxvii_pt::bnw) * x(i+1,ipl-1,k+1);
		}
	}
}


template<>
void copy_rhs<relax_dir::yz>(const stencil_op<seven_pt> & so,
                             const grid_func & x,
                             const grid_func & b,
                             cdr2::mpi::grid_func & b2,
                             int ipl)
{
	for (auto k : b.range(2)) {
		for (auto j : b.range(1)) {
			b2(j,k) = b(ipl,j,k)
				+ so(ipl,j,k,seven_pt::pw) * x(ipl-1,j,k)
				+ so(ipl+1,j,k,seven_pt::pw) * x(ipl+1,j,k);
		}
	}
}


template<>
void copy_rhs<relax_dir::yz>(const stencil_op<xxvii_pt> & so,
                             const grid_func & x,
                             const grid_func & b,
                             cdr2::mpi::grid_func & b2,
                             int ipl)
{
	for (auto k : b.range(2)) {
		for (auto j : b.range(1)) {
			b2(j,k) = b(ipl,j,k)
				+ so(ipl,j+1,k,xxvii_pt::pnw) * x(ipl-1,j+1,k)
				+ so(ipl,j,k,xxvii_pt::pw) * x(ipl-1,j,k)
				+ so(ipl,j,k,xxvii_pt::psw) * x(ipl-1,j-1,k)
				+ so(ipl,j+1,k,xxvii_pt::bnw) * x(ipl-1,j+1,k-1)
				+ so(ipl,j,k,xxvii_pt::bw) * x(ipl-1,j,k-1)
				+ so(ipl,j,k,xxvii_pt::bsw) * x(ipl-1,j-1,k-1)
				+ so(ipl,j+1,k+1,xxvii_pt::bse) * x(ipl-1,j+1,k+1)
				+ so(ipl,j,k+1,xxvii_pt::be) * x(ipl-1,j,k+1)
				+ so(ipl,j,k+1,xxvii_pt::bne) * x(ipl-1,j-1,k+1)
				+ so(ipl+1,j+1,k,xxvii_pt::psw) * x(ipl+1,j+1,k)
				+ so(ipl+1,j,k,xxvii_pt::pw) * x(ipl+1,j,k)
				+ so(ipl+1,j,k,xxvii_pt::pnw) * x(ipl+1,j-1,k)
				+ so(ipl+1,j+1,k,xxvii_pt::bne) * x(ipl+1,j+1,k-1)
				+ so(ipl+1,j,k,xxvii_pt::be) * x(ipl+1,j,k-1)
				+ so(ipl+1,j,k,xxvii_pt::bse) * x(ipl+1,j-1,k-1)
				+ so(ipl+1,j+1,k+1,xxvii_pt::bsw) * x(ipl+1,j+1,k+1)
				+ so(ipl+1,j,k+1,xxvii_pt::bw) * x(ipl+1,j,k+1)
				+ so(ipl+1,j,k+1,xxvii_pt::bnw) * x(ipl+1,j-1,k+1);
		}
	}
}


template<>
void copy32<relax_dir::xy>(const grid_func & x, cdr2::mpi::grid_func & x2, int ipl)
{
	for (auto j : x.grange(1)) {
		for (auto i : x.grange(0)) {
			x2(i,j) = x(i,j,ipl);
		}
	}
}


template<>
void copy23<relax_dir::xy>(const cdr2::mpi::grid_func & x2, grid_func &x, int ipl)
{
	for (auto j : x.grange(1)) {
		for (auto i : x.grange(0)) {
			x(i,j,ipl) = x2(i,j);
		}
	}
}


template<>
void copy32<relax_dir::xz>(const grid_func & x, cdr2::mpi::grid_func & x2, int ipl)
{
	for (auto k : x.grange(2)) {
		for (auto i : x.grange(0)) {
			x2(i,k) = x(i,ipl,k);
		}
	}
}


template<>
void copy23<relax_dir::xz>(const cdr2::mpi::grid_func & x2, grid_func &x, int ipl)
{
	for (auto k : x.grange(2)) {
		for (auto i : x.grange(0)) {
			x(i,ipl,k) = x2(i,k);
		}
	}
}


template<>
void copy32<relax_dir::yz>(const grid_func & x, cdr2::mpi::grid_func & x2, int ipl)
{
	for (auto k : x.grange(2)) {
		for (auto j : x.grange(1)) {
			x2(j,k) = x(ipl,j,k);
		}
	}
}


template<>
void copy23<relax_dir::yz>(const cdr2::mpi::grid_func & x2, grid_func &x, int ipl)
{
	for (auto k : x.grange(2)) {
		for (auto j : x.grange(1)) {
			x(ipl,j,k) = x2(j,k);
		}
	}
}


template<relax_dir rdir>
std::shared_ptr<grid_topo> planes<rdir>::slice_topo(const grid_topo & topo3)
{
	auto igrd = std::make_shared<std::vector<len_t>>(NBMG_pIGRD);
	auto topo2 = std::make_shared<grid_topo>(igrd, 0, 1);

	topo2->nproc(2) = 1;
	if (rdir == relax_dir::xy) {
		MPI_Comm_split(topo3.comm, topo3.coord(2),
		               topo3.coord(1) * topo3.nproc(0) + topo3.coord(0), &topo2->comm);
		for (auto i : range<std::size_t>(2)) {
			topo2->nproc(i) = topo3.nproc(i);
			topo2->coord(i) = topo3.coord(i);
			topo2->is(i) = topo3.is(i);
			topo2->nlocal(i) = topo3.nlocal(i);
			topo2->nglobal(i) = topo3.nglobal(i);
		}

		auto & halo_service = this->services->template get<halo_exchange>();
		auto & dimx = halo_service.leveldims(0);
		auto & dimy = halo_service.leveldims(1);
		topo2->dimxfine.resize(topo2->nproc(0));
		topo2->dimyfine.resize(topo2->nproc(1));
		for (auto i : range<len_t>(topo2->nproc(0))) {
			topo2->dimxfine[i] = dimx(i, topo3.level());
		}

		for (auto j : range<len_t>(topo2->nproc(1))) {
			topo2->dimyfine[j] = dimy(j, topo3.level());
		}
	} else if (rdir == relax_dir::xz) {
		MPI_Comm_split(topo3.comm, topo3.coord(1),
		               topo3.coord(2) * topo3.nproc(0) + topo3.coord(0), &topo2->comm);
		for (auto i : range<std::size_t>(2)) {
			auto i3 = (i == 0) ? 0 : 2;
			topo2->nproc(i) = topo3.nproc(i3);
			topo2->coord(i) = topo3.coord(i3);
			topo2->is(i) = topo3.is(i3);
			topo2->nlocal(i) = topo3.nlocal(i3);
			topo2->nglobal(i) = topo3.nglobal(i3);
		}

		auto & halo_service = this->services->template get<halo_exchange>();
		auto & dimx = halo_service.leveldims(0);
		auto & dimy = halo_service.leveldims(2);
		topo2->dimxfine.resize(topo2->nproc(0));
		topo2->dimyfine.resize(topo2->nproc(1));
		for (auto i : range<len_t>(topo2->nproc(0))) {
			topo2->dimxfine[i] = dimx(i, topo3.level());
		}

		for (auto j : range<len_t>(topo2->nproc(1))) {
			topo2->dimyfine[j] = dimy(j, topo3.level());
		}
	} else if (rdir == relax_dir::yz) {
		MPI_Comm_split(topo3.comm, topo3.coord(0),
		               topo3.coord(2) * topo3.nproc(1) + topo3.coord(1), &topo2->comm);
		for (auto i : range<std::size_t>(2)) {
			auto i3 = i + 1;
			topo2->nproc(i) = topo3.nproc(i3);
			topo2->coord(i) = topo3.coord(i3);
			topo2->is(i) = topo3.is(i3);
			topo2->nlocal(i) = topo3.nlocal(i3);
			topo2->nglobal(i) = topo3.nglobal(i3);
		}

		auto & halo_service = this->services->template get<halo_exchange>();
		auto & dimx = halo_service.leveldims(1);
		auto & dimy = halo_service.leveldims(2);
		topo2->dimxfine.resize(topo2->nproc(0));
		topo2->dimyfine.resize(topo2->nproc(1));
		for (auto i : range<len_t>(topo2->nproc(0))) {
			topo2->dimxfine[i] = dimx(i, topo3.level());
		}

		for (auto j : range<len_t>(topo2->nproc(1))) {
			topo2->dimyfine[j] = dimy(j, topo3.level());
		}
	} else {
		log::error << "invalid relax_dir for planes" << std::endl;
	}

	return topo2;
}


template std::shared_ptr<grid_topo> planes<relax_dir::xy>::slice_topo(const grid_topo & topo3);
template std::shared_ptr<grid_topo> planes<relax_dir::yz>::slice_topo(const grid_topo & topo3);
template std::shared_ptr<grid_topo> planes<relax_dir::xz>::slice_topo(const grid_topo & topo3);



template<relax_dir rdir>
void plane_util<rdir>::copy_coeff(const stencil_op<seven_pt> & so3,
                                  cdr2::mpi::stencil_op<cdr2::five_pt> & so2, int ipl)
{
	using namespace cdr2;

	if (rdir == relax_dir::xy) {
		for (auto j : so3.grange(1)) {
			for (auto i : so3.grange(0)) {
				so2(i,j,five_pt::c) = so3(i,j,ipl,seven_pt::p);
				so2(i,j,five_pt::w) = so3(i,j,ipl,seven_pt::pw);
				so2(i,j,five_pt::s) = so3(i,j,ipl,seven_pt::ps);
			}
		}
	} else if (rdir == relax_dir::xz) {
		for (auto k : so3.grange(2)) {
			for (auto i : so3.grange(0)) {
				so2(i,k,five_pt::c) = so3(i,ipl,k,seven_pt::p);
				so2(i,k,five_pt::w) = so3(i,ipl,k,seven_pt::pw);
				so2(i,k,five_pt::s) = so3(i,ipl,k,seven_pt::b);
			}
		}
	} else if (rdir == relax_dir::yz) {
		for (auto k : so3.grange(2)) {
			for (auto j : so3.grange(1)) {
				so2(j,k,five_pt::c) = so3(ipl,j,k,seven_pt::p);
				so2(j,k,five_pt::w) = so3(ipl,j,k,seven_pt::ps);
				so2(j,k,five_pt::s) = so3(ipl,j,k,seven_pt::b);
			}
		}
	}
}
template void plane_util<relax_dir::xy>::copy_coeff(const stencil_op<seven_pt> & so3,
                                                    cdr2::mpi::stencil_op<cdr2::five_pt> & so2, int ipl);
template void plane_util<relax_dir::yz>::copy_coeff(const stencil_op<seven_pt> & so3,
                                                    cdr2::mpi::stencil_op<cdr2::five_pt> & so2, int ipl);
template void plane_util<relax_dir::xz>::copy_coeff(const stencil_op<seven_pt> & so3,
                                                    cdr2::mpi::stencil_op<cdr2::five_pt> & so2, int ipl);



template<relax_dir rdir>
void plane_util<rdir>::copy_coeff(const stencil_op<xxvii_pt> & so3,
                                  cdr2::mpi::stencil_op<cdr2::nine_pt> & so2, int ipl)
{
	using namespace cdr2;

	if (rdir == relax_dir::xy) {
		for (auto j : so3.grange(1)) {
			for (auto i : so3.grange(0)) {
				so2(i,j,nine_pt::c) = so3(i,j,ipl,xxvii_pt::p);
				so2(i,j,nine_pt::w) = so3(i,j,ipl,xxvii_pt::pw);
				so2(i,j,nine_pt::s) = so3(i,j,ipl,xxvii_pt::ps);
				so2(i,j,nine_pt::sw) = so3(i,j,ipl,xxvii_pt::psw);
				so2(i,j,nine_pt::nw) = so3(i,j,ipl,xxvii_pt::pnw);
			}
		}
	} else if (rdir == relax_dir::xz) {
		for (auto k : so3.grange(2)) {
			for (auto i : so3.grange(0)) {
				so2(i,k,nine_pt::c) = so3(i,ipl,k,xxvii_pt::p);
				so2(i,k,nine_pt::w) = so3(i,ipl,k,xxvii_pt::pw);
				so2(i,k,nine_pt::s) = so3(i,ipl,k,xxvii_pt::b);
				so2(i,k,nine_pt::sw) = so3(i,ipl,k,xxvii_pt::bw);
				so2(i,k,nine_pt::nw) = so3(i,ipl,k,xxvii_pt::be);
			}
		}
	} else if (rdir == relax_dir::yz) {
		for (auto k : so3.grange(2)) {
			for (auto j : so3.grange(1)) {
				so2(j,k,nine_pt::c) = so3(ipl,j,k,xxvii_pt::p);
				so2(j,k,nine_pt::w) = so3(ipl,j,k,xxvii_pt::ps);
				so2(j,k,nine_pt::s) = so3(ipl,j,k,xxvii_pt::b);
				so2(j,k,nine_pt::sw) = so3(ipl,j,k,xxvii_pt::bs);
				so2(j,k,nine_pt::nw) = so3(ipl,j,k,xxvii_pt::bn);
			}
		}
	}
}

template void plane_util<relax_dir::xy>::copy_coeff(const stencil_op<xxvii_pt> & so3,
                                        cdr2::mpi::stencil_op<cdr2::nine_pt> & so2, int ipl);
template void plane_util<relax_dir::yz>::copy_coeff(const stencil_op<xxvii_pt> & so3,
                                        cdr2::mpi::stencil_op<cdr2::nine_pt> & so2, int ipl);
template void plane_util<relax_dir::xz>::copy_coeff(const stencil_op<xxvii_pt> & so3,
                                        cdr2::mpi::stencil_op<cdr2::nine_pt> & so2, int ipl);


template<relax_dir rdir, class sten3, class sten2>
void relax_planes(const stencil_op<sten3> & so, grid_func & x,
                  const grid_func & b, cycle::Dir cdir,
                  services::halo_exchange<stypes> & halo_service,
                  std::vector<slv2_ptr<sten2>> & planes)
{
	auto & conf = planes[0]->get_config();
	auto log_planes = conf.template get<bool>("log-planes", true);
	int lstart, lend, lstride;
	int kgs = so.grid().is(2);

	if (cdir == cycle::Dir::DOWN) {
		lstart = 1 + ((kgs+1) % 2);
		lstride = (kgs % 2) - ((kgs+1) % 2);
		lend = lstart + 2*lstride;
	} else {
		lstart = 1 + (kgs % 2);
		lstride = ((kgs + 1) % 2) - (kgs % 2);
		lend = lstart + 2*lstride;
	}

	// red-black relaxation
	for (auto ipl_beg : range<int>(lstart, lend, lstride)) {
		for (auto ipl : range<int>(ipl_beg, planes.size() + 1, 2)) {
			auto & x2 = planes[ipl-1]->levels.template get<sten2>(0).x;
			auto & b2 = planes[ipl-1]->levels.template get<sten2>(0).b;

			copy32<rdir>(x, x2, ipl);
			copy_rhs<rdir>(so, x, b, b2, ipl);

			auto tmp = log_begin(log_planes, kgs + ipl - 1, relax_dir_name<rdir>::value, x2.grid().comm);
			planes[ipl-1]->solve(b2, x2);
			log_end(log_planes, tmp);

			copy23<rdir>(x2, x, ipl);
		}
		halo_service.run(x, ortho_dir<rdir>::value);
	}
}
template void relax_planes<relax_dir::xy,
                           xxvii_pt, cdr2::nine_pt>(const stencil_op<xxvii_pt> & so, grid_func & x,
                                                    const grid_func & b, cycle::Dir cdir,
                                                    services::halo_exchange<stypes> & halo_service,
                                                    std::vector<slv2_ptr<cdr2::nine_pt>> & planes);
template void relax_planes<relax_dir::yz,
                           xxvii_pt, cdr2::nine_pt>(const stencil_op<xxvii_pt> & so, grid_func & x,
                                                    const grid_func & b, cycle::Dir cdir,
                                                    services::halo_exchange<stypes> & halo_service,
                                                    std::vector<slv2_ptr<cdr2::nine_pt>> & planes);
template void relax_planes<relax_dir::xz,
                           xxvii_pt, cdr2::nine_pt>(const stencil_op<xxvii_pt> & so, grid_func & x,
                                                    const grid_func & b, cycle::Dir cdir,
                                                    services::halo_exchange<stypes> & halo_service,
                                                    std::vector<slv2_ptr<cdr2::nine_pt>> & planes);
template void relax_planes<relax_dir::xy,
                           seven_pt, cdr2::five_pt>(const stencil_op<seven_pt> & so, grid_func & x,
                                                    const grid_func & b, cycle::Dir cdir,
                                                    services::halo_exchange<stypes> & halo_service,
                                                    std::vector<slv2_ptr<cdr2::five_pt>> & planes);
template void relax_planes<relax_dir::yz,
                           seven_pt, cdr2::five_pt>(const stencil_op<seven_pt> & so, grid_func & x,
                                                    const grid_func & b, cycle::Dir cdir,
                                                    services::halo_exchange<stypes> & halo_service,
                                                    std::vector<slv2_ptr<cdr2::five_pt>> & planes);
template void relax_planes<relax_dir::xz,
                           seven_pt, cdr2::five_pt>(const stencil_op<seven_pt> & so, grid_func & x,
                                                    const grid_func & b, cycle::Dir cdir,
                                                    services::halo_exchange<stypes> & halo_service,
                                                    std::vector<slv2_ptr<cdr2::five_pt>> & planes);


template<relax_dir rdir, class sten3, class sten2>
void relax_planes_agg(const stencil_op<sten3> & so, grid_func & x,
                      const grid_func & b, cycle::Dir cdir,
                      services::halo_exchange<stypes> & halo_service,
                      std::vector<slv2_ptr<sten2>> & planes,
                      std::array<plane_ult<sten2>, 2> & threads)
{
	auto & conf = planes[0]->get_config();
	auto log_planes = conf.template get<bool>("log-planes", false);
	int lstart, lend, lstride;
	int kgs = so.grid().is(2);

	if (cdir == cycle::Dir::DOWN) {
		lstart = 1 + ((kgs+1) % 2);
		lstride = (kgs % 2) - ((kgs+1) % 2);
		lend = lstart + 2*lstride;
	} else {
		lstart = 1 + (kgs % 2);
		lstride = ((kgs + 1) % 2) - (kgs % 2);
		lend = lstart + 2*lstride;
	}

	auto tmp = log_begin(log_planes, 0, relax_dir_name<rdir>::value,
	                     planes[0]->levels.template get<sten2>(0).x.grid().comm);
	// red-black relaxation
	for (auto ipl_beg : range<int>(lstart, lend, lstride)) {
		for (auto ipl : range<int>(ipl_beg, planes.size() + 1, 2)) {
			auto & x2 = planes[ipl-1]->levels.template get<sten2>(0).x;
			auto & b2 = planes[ipl-1]->levels.template get<sten2>(0).b;

			copy32<rdir>(x, x2, ipl);
			copy_rhs<rdir>(so, x, b, b2, ipl);

			threads[(ipl-1) % 2].start((ipl - 1) / 2);
		}

		for (auto ipl : range<int>(ipl_beg, planes.size() + 1, 2)) {
			threads[(ipl-1) % 2].join((ipl - 1) / 2);
			auto & x2 = planes[ipl-1]->levels.template get<sten2>(0).x;
			copy23<rdir>(x2, x, ipl);
		}
		halo_service.run(x, ortho_dir<rdir>::value);
	}

	log_end(log_planes, tmp);
}
template void relax_planes_agg<relax_dir::xy,
                               seven_pt, cdr2::five_pt>(const stencil_op<seven_pt> & so, grid_func & x,
                                                        const grid_func & b, cycle::Dir cdir,
                                                        services::halo_exchange<stypes> & halo_service,
                                                        std::vector<slv2_ptr<cdr2::five_pt>> & planes,
                                                        std::array<plane_ult<cdr2::five_pt>, 2> & threads);
template void relax_planes_agg<relax_dir::yz,
                               seven_pt, cdr2::five_pt>(const stencil_op<seven_pt> & so, grid_func & x,
                                                        const grid_func & b, cycle::Dir cdir,
                                                        services::halo_exchange<stypes> & halo_service,
                                                        std::vector<slv2_ptr<cdr2::five_pt>> & planes,
                                                        std::array<plane_ult<cdr2::five_pt>, 2> & threads);
template void relax_planes_agg<relax_dir::xz,
                               seven_pt, cdr2::five_pt>(const stencil_op<seven_pt> & so, grid_func & x,
                                                        const grid_func & b, cycle::Dir cdir,
                                                        services::halo_exchange<stypes> & halo_service,
                                                        std::vector<slv2_ptr<cdr2::five_pt>> & planes,
                                                        std::array<plane_ult<cdr2::five_pt>, 2> & threads);
template void relax_planes_agg<relax_dir::xy,
                               xxvii_pt, cdr2::nine_pt>(const stencil_op<xxvii_pt> & so, grid_func & x,
                                                        const grid_func & b, cycle::Dir cdir,
                                                        services::halo_exchange<stypes> & halo_service,
                                                        std::vector<slv2_ptr<cdr2::nine_pt>> & planes,
                                                        std::array<plane_ult<cdr2::nine_pt>, 2> & threads);
template void relax_planes_agg<relax_dir::yz,
                               xxvii_pt, cdr2::nine_pt>(const stencil_op<xxvii_pt> & so, grid_func & x,
                                                        const grid_func & b, cycle::Dir cdir,
                                                        services::halo_exchange<stypes> & halo_service,
                                                        std::vector<slv2_ptr<cdr2::nine_pt>> & planes,
                                                        std::array<plane_ult<cdr2::nine_pt>, 2> & threads);
template void relax_planes_agg<relax_dir::xz,
                               xxvii_pt, cdr2::nine_pt>(const stencil_op<xxvii_pt> & so, grid_func & x,
                                                        const grid_func & b, cycle::Dir cdir,
                                                        services::halo_exchange<stypes> & halo_service,
                                                        std::vector<slv2_ptr<cdr2::nine_pt>> & planes,
                                                        std::array<plane_ult<cdr2::nine_pt>, 2> & threads);


template<relax_dir rdir>
template<class sten3, class sten2>
void planes<rdir>::setup_impl(const stencil_op<sten3> & so, std::vector<slv2_ptr<sten2>> & planes,
                              std::array<plane_ult<sten2>, 2> & threads,
                              std::array<plane_team, 2> & teams)
{
	#ifdef PLANE_AGG
	this->aggregate = this->params->plane_agg;
	teams[0].threads = threads[0].get_threads();
	teams[1].threads = threads[1].get_threads();
	#else
	this->aggregate = false;
	#endif
	int nplanes = so.shape(2);
	auto rng = so.range(2);
	if (rdir == relax_dir::xz) {
		rng = so.range(1);
		nplanes = so.shape(1);
	} else if (rdir == relax_dir::yz) {
		rng = so.range(0);
		nplanes = so.shape(0);
	}
	auto kgs = so.grid().is(2);
	auto topo2 = slice_topo(so.grid());
	auto conf2 = this->params->plane_config;
	auto log_planes = conf2->template get<bool>("log-planes", true);
	cdr2::mpi::kman_ptr master_kmans[2];
	{
		auto tmp = log_begin(log_planes, kgs + 1 - 1, relax_dir_name<rdir>::value, topo2->comm);
		master_kmans[0] = master_kman(*conf2, (nplanes / 2) + (nplanes % 2), aggregate, teams[0]);
		log_end(log_planes, tmp);
		tmp = log_begin(log_planes, kgs + 2 - 1, relax_dir_name<rdir>::value, topo2->comm);
		master_kmans[1] = master_kman(*conf2, nplanes / 2, aggregate, teams[1]);
		log_end(log_planes, tmp);
	}
	for (auto ipl : rng) {
		int i = ipl-1;
		cdr2::mpi::kman_ptr kman2;
		auto so2_ptr = std::make_unique<cdr2::mpi::stencil_op<sten2>>(topo2);
		auto & so2 = *so2_ptr;
		plane_util<rdir>::copy_coeff(so, so2, ipl);

		auto tmp = log_begin(log_planes, kgs + ipl - 1, relax_dir_name<rdir>::value, topo2->comm);
		if (i < 2)
			kman2 = master_kmans[i];
		else {
			kman2 = worker_kman(*conf2, (i % 2 == 0) ? (nplanes / 2) + (nplanes % 2) : nplanes / 2,
			                    aggregate, teams[i % 2], i / 2);
		}
		planes.emplace_back(std::make_unique<cdr2::mpi::solver<sten2>>(so2, conf2, kman2));
		log_end(log_planes, tmp);

		planes.back()->give_op(std::move(so2_ptr));
		// setup fine-grid solution and right hand side with contiguous memory across planes
		{
			service_manager<cdr2::mpi::stypes> & sman = kman2->services();
			auto & mpool = sman.get<services::mempool>();
			std::size_t nbytes = topo2->nlocal(0) * topo2->nlocal(1) * sizeof(real_t);
			real_t *xaddr = (real_t*) mpool.addr(services::mempool::sol, nbytes);
			real_t *baddr = (real_t*) mpool.addr(services::mempool::rhs, nbytes);

			planes.back()->levels.template get<sten2>(0).x = cdr2::mpi::grid_func(xaddr, topo2);
			planes.back()->levels.template get<sten2>(0).b = cdr2::mpi::grid_func(baddr, topo2);
		}
		#ifdef PLANE_AGG
		if (aggregate)
			threads[i % 2].add_plane(planes.back().get());
		#endif
	}

	// setup services for solve with ults
	#ifdef PLANE_AGG
	if (aggregate) {
		for (auto ipl : rng) {
			int i = ipl - 1;
			setup_agg_solve(*planes[i]);
			planes[i]->apply_heirs([](cdr2::mpi::solver<cdr2::nine_pt> & child) {
					setup_agg_solve(child);
				});
		}
	}
	#endif
}
template void planes<relax_dir::xy>::setup_impl<xxvii_pt, cdr2::nine_pt>(const stencil_op<xxvii_pt> & so, std::vector<slv2_ptr<cdr2::nine_pt>> & planes,
                                                                std::array<plane_ult<cdr2::nine_pt>, 2> & threads,
                                                                std::array<plane_team, 2> & teams);
template void planes<relax_dir::xz>::setup_impl<xxvii_pt, cdr2::nine_pt>(const stencil_op<xxvii_pt> & so, std::vector<slv2_ptr<cdr2::nine_pt>> & planes,
                                                                std::array<plane_ult<cdr2::nine_pt>, 2> & threads,
                                                                std::array<plane_team, 2> & teams);
template void planes<relax_dir::yz>::setup_impl<xxvii_pt, cdr2::nine_pt>(const stencil_op<xxvii_pt> & so, std::vector<slv2_ptr<cdr2::nine_pt>> & planes,
                                                                std::array<plane_ult<cdr2::nine_pt>, 2> & threads,
                                                                std::array<plane_team, 2> & teams);
template void planes<relax_dir::xy>::setup_impl<seven_pt, cdr2::five_pt>(const stencil_op<seven_pt> & so, std::vector<slv2_ptr<cdr2::five_pt>> & planes,
                                                                std::array<plane_ult<cdr2::five_pt>, 2> & threads,
                                                                std::array<plane_team, 2> & teams);
template void planes<relax_dir::xz>::setup_impl<seven_pt, cdr2::five_pt>(const stencil_op<seven_pt> & so, std::vector<slv2_ptr<cdr2::five_pt>> & planes,
                                                                std::array<plane_ult<cdr2::five_pt>, 2> & threads,
                                                                std::array<plane_team, 2> & teams);
template void planes<relax_dir::yz>::setup_impl<seven_pt, cdr2::five_pt>(const stencil_op<seven_pt> & so, std::vector<slv2_ptr<cdr2::five_pt>> & planes,
                                                                std::array<plane_ult<cdr2::five_pt>, 2> & threads,
                                                                std::array<plane_team, 2> & teams);



}}}
