#ifndef CEDAR_3D_MPI_RELAX_PLANES_H
#define CEDAR_3D_MPI_RELAX_PLANES_H

#include <tuple>

#include <cedar/types.h>
#include <cedar/kernels/plane_relax.h>
#include <cedar/3d/mpi/types.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/3d/mpi/kernel_manager.h>
#include <cedar/3d/mpi/plane_ult.h>
#include <cedar/3d/mpi/plane_team.h>

namespace cedar { namespace cdr3 { namespace mpi {

std::tuple<int,MPI_Comm> log_begin(bool log_planes, int ipl, const std::string & suff, MPI_Comm comm);
void log_end(bool log_planes, std::tuple<int,MPI_Comm> saved);

cdr2::mpi::kman_ptr master_kman(config & conf, int nplanes, bool aggregate, plane_team & team);
cdr2::mpi::kman_ptr worker_kman(config & conf, int nplanes, bool aggregate, plane_team & team, int worker_id);

template<class fsten>
void setup_agg_solve(cdr2::mpi::solver<fsten> & slv)
{
	service_manager<cdr2::mpi::stypes> & sman = slv.get_kernels()->services();
	sman.set<services::halo_exchange<cdr2::mpi::stypes>>("plane");
	sman.set<services::message_passing>("plane");
	slv.setup_halo();
}


template<class fsten>
using slv2_ptr = std::unique_ptr<cdr2::mpi::solver<fsten>>;

template<relax_dir rdir, class sten3>
void copy_rhs(const stencil_op<sten3> & so,
              const grid_func & x,
              const grid_func & b,
              cdr2::mpi::grid_func & b2,
              int ipl);

template<relax_dir rdir>
void copy32(const grid_func & x, cdr2::mpi::grid_func & x2, int ipl);

template<relax_dir rdir>
void copy23(const cdr2::mpi::grid_func & x2, grid_func &x, int ipl);


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
		halo_service.run(x, 4);
	}
}


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
		halo_service.run(x, 4);
	}

	log_end(log_planes, tmp);
}


template<relax_dir rdir, class sten3, class sten2>
void copy_coeff(const stencil_op<sten3> & so3,
                cdr2::mpi::stencil_op<sten2> & so2);

template<relax_dir rdir>
void copy_coeff(const stencil_op<seven_pt> & so3,
                cdr2::mpi::stencil_op<cdr2::five_pt> & so2)
{
	using namespace cdr2;

	if (rdir == relax_dir::xy) {
		for (auto k : so3.range(2)) {
			for (auto j : so3.grange(1)) {
				for (auto i : so3.grange(0)) {
					so2(i,j,five_pt::c) = so3(i,j,k,seven_pt::p);
					so2(i,j,five_pt::w) = so3(i,j,k,seven_pt::pw);
					so2(i,j,five_pt::s) = so3(i,j,k,seven_pt::ps);
				}
			}
		}
	} else if (rdir == relax_dir::xz) {
		for (auto j : so3.range(1)) {
			for (auto k : so3.grange(2)) {
				for (auto i : so3.grange(0)) {
					so2(i,k,five_pt::c) = so3(i,j,k,seven_pt::p);
					so2(i,k,five_pt::w) = so3(i,j,k,seven_pt::pw);
					so2(i,k,five_pt::s) = so3(i,j,k,seven_pt::b);
				}
			}
		}
	} else if (rdir == relax_dir::yz) {
		for (auto i : so3.range(0)) {
			for (auto k : so3.grange(2)) {
				for (auto j : so3.grange(1)) {
					so2(j,k,five_pt::c) = so3(i,j,k,seven_pt::p);
					so2(j,k,five_pt::w) = so3(i,j,k,seven_pt::ps);
					so2(j,k,five_pt::s) = so3(i,j,k,seven_pt::b);
				}
			}
		}
	}
}

template<relax_dir rdir>
void copy_coeff(const stencil_op<xxvii_pt> & so3,
                cdr2::mpi::stencil_op<cdr2::nine_pt> & so2)
{
	using namespace cdr2;

	if (rdir == relax_dir::xy) {
		for (auto k : so3.range(2)) {
			for (auto j : so3.grange(1)) {
				for (auto i : so3.grange(0)) {
					so2(i,j,nine_pt::c) = so3(i,j,k,xxvii_pt::p);
					so2(i,j,nine_pt::w) = so3(i,j,k,xxvii_pt::pw);
					so2(i,j,nine_pt::s) = so3(i,j,k,xxvii_pt::ps);
					so2(i,j,nine_pt::sw) = so3(i,j,k,xxvii_pt::psw);
					so2(i,j,nine_pt::nw) = so3(i,j,k,xxvii_pt::pnw);
				}
			}
		}
	} else if (rdir == relax_dir::xz) {
		for (auto j : so3.range(1)) {
			for (auto k : so3.grange(2)) {
				for (auto i : so3.grange(0)) {
					so2(i,k,nine_pt::c) = so3(i,j,k,xxvii_pt::p);
					so2(i,k,nine_pt::w) = so3(i,j,k,xxvii_pt::pw);
					so2(i,k,nine_pt::s) = so3(i,j,k,xxvii_pt::b);
					so2(i,k,nine_pt::sw) = so3(i,j,k,xxvii_pt::bw);
					so2(i,k,nine_pt::nw) = so3(i,j,k,xxvii_pt::be);
				}
			}
		}
	} else if (rdir == relax_dir::yz) {
		for (auto i : so3.range(0)) {
			for (auto k : so3.grange(2)) {
				for (auto j : so3.grange(1)) {
					so2(j,k,nine_pt::c) = so3(i,j,k,xxvii_pt::p);
					so2(j,k,nine_pt::w) = so3(i,j,k,xxvii_pt::ps);
					so2(j,k,nine_pt::s) = so3(i,j,k,xxvii_pt::b);
					so2(j,k,nine_pt::sw) = so3(i,j,k,xxvii_pt::bs);
					so2(j,k,nine_pt::nw) = so3(i,j,k,xxvii_pt::bn);
				}
			}
		}
	}
}

template<relax_dir rdir>
class planes : public kernels::plane_relax<stypes, rdir>
{
public:
	planes() : aggregate(false) {}

	void setup(const stencil_op<seven_pt> & so) override
	{
		this->setup_impl(so, fine_planes, fine_threads);
	}
	void setup(const stencil_op<xxvii_pt> & so) override
	{
		level_threads.emplace_back();
		level_planes.emplace_back();
		level_map[so.shape(0)] = level_planes.size() - 1;
		this->setup_impl(so, level_planes.back(), level_threads.back());
	}
	void run(const stencil_op<seven_pt> & so, grid_func & x,
	         const grid_func & b, cycle::Dir cycle_dir) override { this->run_impl(so, x, b, cycle_dir); }
	void run(const stencil_op<xxvii_pt> & so, grid_func & x,
	         const grid_func & b, cycle::Dir cycle_dir) override { this->run_impl(so, x, b, cycle_dir); }

	template<class sten3, class sten2>
	void setup_impl(const stencil_op<sten3> & so, std::vector<slv2_ptr<sten2>> & planes,
	                std::array<plane_ult<sten2>, 2> & threads)
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
			copy_coeff<rdir>(so, so2);

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

	template<class sten>
	void run_impl(const stencil_op<sten> & so, grid_func & x,
	              const grid_func & b, cycle::Dir cycle_dir)
	{
		auto & halo_service = this->services->template get<halo_exchange>();
		if (std::is_same<sten, seven_pt>::value) {
			#ifdef PLANE_AGG
			if (aggregate)
				relax_planes_agg<rdir>(so, x, b, cycle_dir, halo_service, fine_planes, fine_threads);
			else
				relax_planes<rdir>(so, x, b, cycle_dir, halo_service, fine_planes);
			#else
			relax_planes<rdir>(so, x, b, cycle_dir, halo_service, fine_planes);
			#endif
		} else {
			len_t lsize = so.shape(0);
			#ifdef PLANE_AGG
			if (aggregate)
				relax_planes_agg<rdir>(so, x, b, cycle_dir, halo_service, level_planes[level_map[lsize]], level_threads[level_map[lsize]]);
			else
				relax_planes<rdir>(so, x, b, cycle_dir, halo_service, level_planes[level_map[lsize]]);
			#else
			relax_planes<rdir>(so, x, b, cycle_dir, halo_service, level_planes[level_map[lsize]]);
			#endif
		}
	}

protected:
	std::vector<std::vector<slv2_ptr<cdr2::nine_pt>>> level_planes;
	std::vector<slv2_ptr<cdr2::five_pt>> fine_planes;
	std::vector<std::array<plane_ult<cdr2::nine_pt>, 2>> level_threads;
	std::array<plane_ult<cdr2::five_pt>, 2> fine_threads;
	std::map<len_t, std::size_t> level_map;
	std::array<plane_team, 2> teams;
	bool aggregate;

	std::shared_ptr<grid_topo> slice_topo(const grid_topo & topo3)
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
};

}}}

#endif
