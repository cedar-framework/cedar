#ifndef CEDAR_2D_SOLVER_MPI_CEDAR_H
#define CEDAR_2D_SOLVER_MPI_CEDAR_H

#include <algorithm>
#include <memory>
#include <mpi.h>
#include <array>

#include <cedar/multilevel.h>
#include <cedar/level.h>
#include <cedar/2d/level_container.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/solver.h>
#include <cedar/2d/mpi/types.h>
#include <cedar/2d/mpi/kernel_manager.h>
#include <cedar/perf/predict.h>
#include <cedar/2d/mpi/redist_solver.h>
#include <cedar/2d/redist/multilevel_wrapper.h>
#include <cedar/2d/redist/cholesky_solver.h>
#include <cedar/2d/redist/setup_redist.h>

namespace cedar { namespace cdr2 { namespace mpi {

template<class sten>
	struct level2mpi : public level<sten, stypes>
{
	using parent = level<sten, stypes>;
	level2mpi(services::mempool & mpool, topo_ptr topo) : parent::level(topo)
	{
		std::size_t nbytes = topo->nlocal(0)*topo->nlocal(1)*sizeof(real_t);
		real_t *xaddr   = (real_t*) mpool.addr(mempool::sol, nbytes);
		real_t *resaddr = (real_t*) mpool.addr(mempool::res, nbytes);
		real_t *baddr   = (real_t*) mpool.addr(mempool::rhs, nbytes);

		this->x = grid_func(xaddr, topo);
		this->res = grid_func(resaddr, topo);
		this->b = grid_func(baddr, topo);
		this->SOR = {{relax_stencil(topo->nlocal(0)-2, topo->nlocal(1)-2),
		              relax_stencil(topo->nlocal(0)-2, topo->nlocal(1)-2)}};
		this->R.associate(&this->P);
	}


	level2mpi(services::mempool & mpool, stencil_op<sten> & A) : parent::level(A)
	{
		topo_ptr topo = A.grid_ptr();
		std::size_t nbytes = topo->nlocal(0)*topo->nlocal(1)*sizeof(real_t);
		this->res = mpi::grid_func((real_t*) mpool.addr(mempool::res, nbytes),
		                           topo);
		this->SOR = {{relax_stencil(A.shape(0), A.shape(1)),
		              relax_stencil(A.shape(0), A.shape(1))}};
	}
};


template<class fsten>
	class solver: public multilevel<exec_mode::mpi, level_container<level2mpi, fsten>,
	fsten, cdr2::mpi::solver<fsten>>
{
public:
	using parent = multilevel<exec_mode::mpi, level_container<level2mpi, fsten>,
	                          fsten, cdr2::mpi::solver<fsten>>;
	using parent::kman;
solver(mpi::stencil_op<fsten> & fop) : parent::multilevel(fop), comm(fop.grid().comm)
	{
		this->kman = build_kernel_manager(*this->conf);
		parent::setup(fop);
	}


	solver(mpi::stencil_op<fsten> & fop,
	       std::shared_ptr<config> conf) : parent::multilevel(fop, conf), comm(fop.grid().comm)
	{
		this->kman = build_kernel_manager(*this->conf);
		parent::setup(fop);
	}


	solver(mpi::stencil_op<fsten> & fop,
	       std::shared_ptr<config> conf,
	       kman_ptr kman) : parent::multilevel(fop, conf), comm(fop.grid().comm)
	{
		this->kman = kman;
		parent::setup(fop);
	}


	~solver() {
	}


	std::size_t compute_num_levels(mpi::stencil_op<fsten> & fop)
	{
		int ng;

		kman->template run<setup_nog>(fop.grid(), this->settings.min_coarse, &ng);

		return ng;
	}


	virtual cdr2::mpi::grid_func solve(const cdr2::mpi::grid_func & b) override
	{
		auto & bd = const_cast<cdr2::mpi::grid_func&>(b);
		auto & halo_service = kman->services().template get<halo_exchange>();
		halo_service.run(bd);
		return parent::solve(b);
	}


	virtual void solve(const cdr2::mpi::grid_func & b, cdr2::mpi::grid_func & x) override
	{
		auto & bd = const_cast<cdr2::mpi::grid_func&>(b);
		auto & halo_service = kman->services().template get<halo_exchange>();
		halo_service.run(bd);
		return parent::solve(b, x);
	}

	virtual void setup_cg_solve() override
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
			if ((choice[0] != 1) and (choice[1] != 1)) {
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


	grid_topo & get_grid(std::size_t i)
	{
		if (i == 0) {
			auto & fop = this->levels.template get<fsten>(i).A;
			return fop.grid();
		} else {
			auto & sop = this->levels.get(i).A;
			return sop.grid();
		}
	}


	void setup_space(std::size_t nlevels)
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
	}

	void setup_halo()
	{
		auto & sop = this->levels.template get<fsten>(0).A;

		std::vector<topo_ptr> topos;
		topos.push_back(sop.grid_ptr());

		for (std::size_t i = 1; i < this->nlevels(); i++)
			topos.push_back(this->levels.get(i).A.grid_ptr());

		auto & halo_service = kman->services().template get<halo_exchange>();
		halo_service.setup(topos);
	}

	void give_op(std::unique_ptr<stencil_op<fsten>> fop) {fop_ref = std::move(fop);}


	void apply_heirs(std::function<void(solver<nine_pt> &)> fun)
	{
		if (heir) {
			auto & slv = heir->get_inner().get_inner();
			fun(slv);
			slv.apply_heirs(fun);
		}
	}


	MPI_Comm comm;


protected:
	std::unique_ptr<stencil_op<fsten>> fop_ref;
	std::shared_ptr<redist_solver<multilevel_wrapper<solver<nine_pt>>>> heir;
};

}}}

#endif
