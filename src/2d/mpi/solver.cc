#include <cedar/2d/mpi/solver.h>
#include <ftl/Device.hpp>

namespace cedar { namespace cdr2 { namespace mpi {

template<class fsten>
std::size_t solver<fsten>::compute_num_levels(stencil_op<fsten> & fop)
{
	int ng;

	kman->template run<setup_nog>(fop.grid(), this->settings.min_coarse, &ng);

	return ng;
}
template std::size_t solver<five_pt>::compute_num_levels(stencil_op<five_pt> & fop);
template std::size_t solver<nine_pt>::compute_num_levels(stencil_op<nine_pt> & fop);


template<class fsten>
void solver<fsten>::solve(const cdr2::mpi::grid_func & b, cdr2::mpi::grid_func & x)
{
	auto & bd = const_cast<cdr2::mpi::grid_func&>(b);
	auto & halo_service = kman->services().template get<halo_exchange>();
	halo_service.run(bd);
	return parent::solve(b, x);
}
template void solver<five_pt>::solve(const cdr2::mpi::grid_func & b, cdr2::mpi::grid_func & x);
template void solver<nine_pt>::solve(const cdr2::mpi::grid_func & b, cdr2::mpi::grid_func & x);


template<class fsten>
grid_func solver<fsten>::solve(const grid_func & b)
{
	auto & bd = const_cast<cdr2::mpi::grid_func&>(b);
	auto & halo_service = kman->services().template get<halo_exchange>();
	halo_service.run(bd);
	return parent::solve(b);
}
template grid_func solver<five_pt>::solve(const grid_func & b);
template grid_func solver<nine_pt>::solve(const grid_func & b);


template<class fsten>
void solver<fsten>::setup_cg_solve()
{
	auto params = build_kernel_params(*(this->conf));
	auto & cop = this->levels.get(this->levels.size() - 1).A;

	auto cg_conf = this->settings.coarse_config;

	if (this->settings.coarse_solver == ml_settings::cg_type::redist) {
		auto & fgrid = this->levels.template get<fsten>(0).A.grid();
		auto choice = choose_redist<2>(this->settings.rsettings,
		                               std::array<int,2>({{fgrid.nproc(0), fgrid.nproc(1)}}),
		                               std::array<len_t,2>({{fgrid.nglobal(0), fgrid.nglobal(1)}}));
		MPI_Bcast(choice.data(), 2, MPI_INT, 0, fgrid.comm);
		if (not (choice[0]*choice[1] == 1)) {
			log::status << "Redistributing to " << choice[0] << " x " << choice[1] << " cores" << std::endl;
			using inner_solver = multilevel_wrapper<mpi::solver<nine_pt>>;
			this->heir = create_redist_ptr<inner_solver>(kman, cop, cg_conf, choice);
			this->coarse_solver = wrap_redist_ptr<inner_solver>(*this->conf, kman, this->heir);
			return;
		}
	}

	std::array<int, 2> choice{{1,1}};

	log::status << "Redistributing to " << choice[0] << " x " << choice[1] << " cores" << std::endl;
	if (this->settings.coarse_solver == ml_settings::cg_type::lu) {
		this->coarse_solver = create_redist_solver<cholesky_solver>(kman,
		                                                            *this->conf,
		                                                            cop,
		                                                            cg_conf,
		                                                            choice);
	} else {
		using inner_solver = multilevel_wrapper<cdr2::solver<nine_pt>>;
		this->coarse_solver = create_redist_solver<inner_solver>(kman,
		                                                         *this->conf,
		                                                         cop,
		                                                         cg_conf,
		                                                         choice);
	}
}
template void solver<five_pt>::setup_cg_solve();
template void solver<nine_pt>::setup_cg_solve();


template<class fsten>
void solver<fsten>::setup_space(std::size_t nlevels)
{
	service_manager<stypes> & sman = kman->services();
	this->levels.init(sman.get<mempool>(),
	                  nlevels);
	for (auto i : range<std::size_t>(nlevels-1)) {

		auto & fgrid = this->get_grid(i);
		if (i == 0)
			fgrid.grow(nlevels);

		int kc = nlevels - i - 2;

		auto cgrid = std::make_shared<grid_topo>(fgrid.get_igrd(), kc, nlevels);
		cgrid->comm = fgrid.comm;

		len_t NLxg = fgrid.nlocal(0) - 2;
		len_t NLyg = fgrid.nlocal(1) - 2;
		len_t NGxg = (fgrid.nglobal(0) - 1) / 2 + 2;
		len_t NGyg = (fgrid.nglobal(1) - 1) / 2 + 2;

		cgrid->nglobal(0) = NGxg;
		cgrid->nglobal(1) = NGyg;

		if ((fgrid.is(0) % 2) == 1) {
			cgrid->is(0) = (fgrid.is(0) + 1) / 2;
			NLxg = (NLxg + 1) / 2;
		} else {
			cgrid->is(0) = fgrid.is(0)/2 + 1;
			if (NLxg % 2 == 1) NLxg = (NLxg-1)/2;
			else NLxg = (NLxg+1)/2;
		}


		if (fgrid.is(1) % 2 == 1) {
			cgrid->is(1) = (fgrid.is(1)+1) / 2;
			NLyg = (NLyg+1) / 2;
		} else {
			cgrid->is(1) = fgrid.is(1) / 2 + 1;
			if (NLyg % 2 == 1) NLyg = (NLyg - 1) / 2;
			else NLyg = (NLyg+1)/2;
		}

		cgrid->nlocal(0) = NLxg + 2;
		cgrid->nlocal(1) = NLyg + 2;

		cgrid->nproc(0) = fgrid.nproc(0);
		cgrid->nproc(1) = fgrid.nproc(1);
		cgrid->coord(0) = fgrid.coord(0);
		cgrid->coord(1) = fgrid.coord(1);

		this->levels.add(sman.get<mempool>(),
		                 cgrid);
	}
	setup_halo();
	{
		auto & sop = this->levels.template get<fsten>(0).A;
		auto & halo_service = kman->services().template get<halo_exchange>();
		halo_service.run(sop);
	}

        /* Do GPU setup if it's enabled */
        if (this->settings.use_gpu) {
            ftl::device::autodetect();
            log::status << "Running on " << ftl::device::get_name() << std::endl;

            /* Load kernels */
            if (!ftl::load_kernels(true)) {
                throw std::runtime_error("FTL layer: Failed to load kernels.");
            }
            log::status << "Loaded and compiled all GPU kernels" << std::endl;
            log::status << "Allocating GPU memory" << std::endl;

            for (std::size_t i = 0; i < this->nlevels(); i++) {
                if (i == 0) {
                    auto& level = this->levels.template get<fsten>(i);
                    level.A.ensure_gpu();
                    level.P.ensure_gpu();
                    level.x.ensure_gpu();
                    level.res.ensure_gpu();
                    level.b.ensure_gpu();
                } else {
                    auto& level = this->levels.get(i);
                    level.A.ensure_gpu();
                    level.P.ensure_gpu();
                    level.x.ensure_gpu();
                    level.res.ensure_gpu();
                    level.b.ensure_gpu();
                }
            }
        }
}
template void solver<five_pt>::setup_space(std::size_t nlevels);
template void solver<nine_pt>::setup_space(std::size_t nlevels);


template<class fsten>
void solver<fsten>::setup_halo()
{
	auto & sop = this->levels.template get<fsten>(0).A;

	std::vector<topo_ptr> topos;
	topos.push_back(sop.grid_ptr());

	for (std::size_t i = 1; i < this->nlevels(); i++)
		topos.push_back(this->levels.get(i).A.grid_ptr());

	service_manager<stypes> & sman = kman->services();
	auto & halo_service = sman.get<halo_exchange>();
	halo_service.setup(topos);
}
template void solver<five_pt>::setup_halo();
template void solver<nine_pt>::setup_halo();


template<class fsten>
void solver<fsten>::apply_heirs(std::function<void(solver<nine_pt> &)> fun)
{
	if (heir and (heir->isactive())) {
		auto & slv = heir->get_inner().get_inner();
		fun(slv);
		slv.apply_heirs(fun);
	}
}
template void solver<five_pt>::apply_heirs(std::function<void(solver<nine_pt> &)> fun);
template void solver<nine_pt>::apply_heirs(std::function<void(solver<nine_pt> &)> fun);

}}}
