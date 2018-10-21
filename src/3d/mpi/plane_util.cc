#include <cedar/3d/mpi/plane_mempool.h>
#include <cedar/3d/mpi/plane_mpi.h>
#include <cedar/3d/mpi/plane_exchange.h>
#include <cedar/3d/mpi/plane_util.h>


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
	timer_down();
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

			timer_begin("plane-copy");
			copy32<rdir>(x, x2, ipl);
			copy_rhs<rdir>(so, x, b, b2, ipl);
			timer_end("plane-copy");

			auto tmp = log_begin(log_planes, kgs + ipl - 1, relax_dir_name<rdir>::value, x2.grid().comm);
			planes[ipl-1]->solve(b2, x2);
			log_end(log_planes, tmp);

			timer_begin("plane-copy");
			copy23<rdir>(x2, x, ipl);
			timer_end("plane-copy");
		}
		timer_begin("halo");
		halo_service.run(x, ortho_dir<rdir>::value);
		timer_end("halo");
	}
	timer_up();
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
	timer_down();
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

			timer_begin("plane-copy");
			copy32<rdir>(x, x2, ipl);
			copy_rhs<rdir>(so, x, b, b2, ipl);
			timer_end("plane-copy");

			threads[(ipl-1) % 2].start((ipl - 1) / 2);
		}

		for (auto ipl : range<int>(ipl_beg, planes.size() + 1, 2)) {
			threads[(ipl-1) % 2].join((ipl - 1) / 2);
			auto & x2 = planes[ipl-1]->levels.template get<sten2>(0).x;
			timer_begin("plane-copy");
			copy23<rdir>(x2, x, ipl);
			timer_end("plane-copy");
		}
		timer_begin("halo");
		halo_service.run(x, ortho_dir<rdir>::value);
		timer_end("halo");
	}

	log_end(log_planes, tmp);
	timer_up();
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
}}}
