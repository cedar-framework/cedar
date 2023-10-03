#include <cedar/3d/mpi/plane_mempool.h>
#include <cedar/3d/mpi/plane_mpi.h>
#include <cedar/3d/mpi/relax_planes.h>


namespace cedar { namespace cdr3 { namespace mpi {

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
	(void)threads; // avoid warning
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
	auto log_planes = conf2->template get<bool>("log-planes", false);
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

		timer_pause();
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
		timer_play();
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
