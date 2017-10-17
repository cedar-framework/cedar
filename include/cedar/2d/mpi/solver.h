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
#include <cedar/2d/relax_stencil.h>
#include <cedar/2d/inter/mpi/prolong_op.h>
#include <cedar/2d/inter/mpi/restrict_op.h>
#include <cedar/2d/kernel/mpi/registry.h>
#include <cedar/2d/mpi/msg_exchanger.h>
#include <cedar/2d/solver.h>
#include <cedar/perf/predict.h>
#include <cedar/2d/mpi/redist_solver.h>
#include <cedar/2d/redist/multilevel_wrapper.h>

namespace cedar { namespace cdr2 { namespace mpi {

template<class sten>
	struct level2mpi : public level<sten, stypes>
{
	using parent = level<sten, stypes>;
level2mpi(topo_ptr topo) : parent::level(topo)
	{
		this->SOR = {{relax_stencil(topo->nlocal(0)-2, topo->nlocal(1)-2),
		              relax_stencil(topo->nlocal(0)-2, topo->nlocal(1)-2)}};
		this->R.associate(&this->P);
	}


level2mpi(stencil_op<sten> & A) : parent::level(A)
	{
		this->res = mpi::grid_func(A.grid_ptr());
		this->SOR = {{relax_stencil(A.shape(0), A.shape(1)),
		              relax_stencil(A.shape(0), A.shape(1))}};
	}
};


template<class fsten>
	class solver: public multilevel<level_container<level2mpi, fsten>,
	typename kernel::mpi::registry::parent, fsten, cdr2::mpi::solver<fsten>>
{
public:
	using parent = multilevel<level_container<level2mpi, fsten>,
		typename kernel::mpi::registry::parent, fsten, cdr2::mpi::solver<fsten>>;
solver(mpi::stencil_op<fsten> & fop) : parent::multilevel(fop), comm(fop.grid().comm)
	{
		this->kreg = std::make_shared<kernel::mpi::registry>(*(this->conf));
		parent::setup(fop);
	}


	solver(mpi::stencil_op<fsten> & fop,
	       std::shared_ptr<config::reader> conf) : parent::multilevel(fop, conf), comm(fop.grid().comm)
	{
		this->kreg = std::make_shared<kernel::mpi::registry>(*(this->conf));
		parent::setup(fop);
	}


	~solver() {
		// TODO: what??
		if (cg_solver_lu) this->bbd = new real_t[1];
	}


	std::size_t compute_num_levels(mpi::stencil_op<fsten> & fop)
	{
		int ng;
		auto min_coarse = this->conf->template get<len_t>("solver.min_coarse", 3);
		auto kernels = this->kernel_registry();

		kernels->setup_nog(fop.grid(), min_coarse, &ng);

		return ng;
	}


	virtual cdr2::mpi::grid_func solve(const cdr2::mpi::grid_func & b) override
	{
		auto kernels = this->kernel_registry();
		auto & bd = const_cast<cdr2::mpi::grid_func&>(b);
		kernels->halo_exchange(bd);
		return parent::solve(b);
	}


	virtual void solve(const cdr2::mpi::grid_func & b, cdr2::mpi::grid_func & x) override
	{
		auto kernels = this->kernel_registry();
		auto & bd = const_cast<cdr2::mpi::grid_func&>(b);
		kernels->halo_exchange(bd);
		return parent::solve(b, x);
	}

	virtual void setup_cg_solve() override
	{
		auto params = build_kernel_params(*(this->conf));
		auto kernels = this->kernel_registry();
		auto & cop = this->levels.get(this->levels.size() - 1).A;
		std::string cg_solver_str = this->conf->template get<std::string>("solver.cg-solver", "LU");

		auto cg_conf = this->conf->getconf("cg-config");
		if (!cg_conf)
			cg_conf = this->conf;

		if (cg_solver_str != "LU") {
			auto & fgrid = this->levels.template get<fsten>(0).A.grid();
			auto choice = choose_redist<2>(*(this->conf),
			                               std::array<int,2>({fgrid.nproc(0), fgrid.nproc(1)}),
			                               std::array<len_t,2>({fgrid.nglobal(0), fgrid.nglobal(1)}));
			MPI_Bcast(choice.data(), 2, MPI_INT, 0, fgrid.comm);
			if ((choice[0] != 1) and (choice[1] != 1)) {
				log::status << "Redistributing to " << choice[0] << " x " << choice[1] << " cores" << std::endl;
				auto cg_bmg = std::make_shared<mpi::redist_solver<multilevel_wrapper<mpi::solver<nine_pt>>>>(cop,
					                       kernels->get_halo_exchanger().get(),
					                       cg_conf,
					                       choice);
				this->coarse_solver = [&,cg_bmg,kernels,params](mpi::grid_func &x, const mpi::grid_func &b) {
					cg_bmg->solve(x, b);
					if (params->per_mask())
						kernels->halo_exchange(x);
				};
				return;
			}
		}

		std::array<int, 2> choice{{1,1}};

		log::status << "Redistributing to " << choice[0] << " x " << choice[1] << " cores" << std::endl;
		auto cg_bmg = std::make_shared<
			mpi::redist_solver<
			multilevel_wrapper<
			cdr2::solver<nine_pt>
			>>>(cop,
			    kernels->get_halo_exchanger().get(),
			    cg_conf,
			    choice);
		this->coarse_solver = [&,cg_bmg,kernels,params](mpi::grid_func &x, const mpi::grid_func &b) {
			cg_bmg->solve(x, b);
			if (params->per_mask())
				kernels->halo_exchange(x);
		};
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
		this->levels.init(nlevels);
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

			this->levels.add(cgrid);
		}
		setup_halo();
	}

	void setup_halo()
	{
		auto & sop = this->levels.template get<fsten>(0).A;

		std::vector<topo_ptr> topos;
		topos.push_back(sop.grid_ptr());

		for (std::size_t i = 1; i < this->nlevels(); i++)
			topos.push_back(this->levels.get(i).A.grid_ptr());

		this->kreg->halo_setup(topos);
		this->kreg->halo_stencil_exchange(sop);
	}

	MPI_Comm comm;

private:
	bool cg_solver_lu;
};

}}}

#endif
