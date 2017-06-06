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
#include <cedar/2d/mpi/halo.h>
#include <cedar/2d/solver.h>
#include <cedar/perf/predict.h>

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
		kernels->halo_exchange(bd, halo_ctx);
		return parent::solve(b);
	}


	virtual void solve(const cdr2::mpi::grid_func & b, cdr2::mpi::grid_func & x) override
	{
		auto kernels = this->kernel_registry();
		auto & bd = const_cast<cdr2::mpi::grid_func&>(b);
		kernels->halo_exchange(bd, halo_ctx);
		return parent::solve(b, x);
	}

	virtual void setup_cg_solve() override
	{
		auto params = build_kernel_params(*(this->conf));
		auto kernels = this->kernel_registry();
		auto & cop = this->levels.get(this->levels.size() - 1).A;
		std::string cg_solver_str = this->conf->template get<std::string>("solver.cg-solver", "LU");

		if (cg_solver_str == "LU")
			cg_solver_lu = true;
		else
			cg_solver_lu = false;

		if (cg_solver_lu) {
			auto & coarse_topo = cop.grid();
			auto nxc = coarse_topo.nglobal(0);
			auto nyc = coarse_topo.nglobal(1);
			this->ABD = mpi::grid_func(nxc+2, nxc*nyc);
			this->bbd = new real_t[this->ABD.len(1)];
			parent::setup_cg_solve();
		} else {
			auto cg_conf = this->conf->getconf("cg-config");
			if (!cg_conf)
				cg_conf = this->conf;

			if (cg_solver_str == "redist") {
				// TODO: should this be the coarse grid!
				auto & fgrid = this->levels.template get(0).A.grid();
				auto choice = choose_redist<2>(*(this->conf),
				                               std::array<int,2>({fgrid.nproc(0), fgrid.nproc(1)}),
				                               std::array<len_t,2>({fgrid.nglobal(0), fgrid.nglobal(1)}));
				MPI_Bcast(choice.data(), 2, MPI_INT, 0, fgrid.comm);
				if ((choice[0] == 1) or (choice[1] == 1)) {
					cg_solver_str = "boxmg";
				} else {
					log::status << "Redistributing to " << choice[0] << " x " << choice[1] << " cores" << std::endl;
					std::shared_ptr<mpi::redist_solver> cg_bmg;
					std::vector<int> nblocks(choice.begin(), choice.end());
					kernels->setup_cg_redist(cop, cg_conf, &cg_bmg, nblocks);
					this->coarse_solver = [&,cg_bmg,kernels,params](mpi::grid_func &x, const mpi::grid_func &b) {
						kernels->solve_cg_redist(*cg_bmg, x, b);
						if (params->per_mask())
							kernels->halo_exchange(x);
					};
				}
			}

			if (cg_solver_str == "boxmg") {
				std::shared_ptr<cdr2::solver<nine_pt>> cg_bmg;
				kernels->setup_cg_boxmg(cop, cg_conf, &cg_bmg);
				this->coarse_solver = [&,cg_bmg,kernels,params](mpi::grid_func &x, const mpi::grid_func &b) {
					kernels->solve_cg_boxmg(*cg_bmg, x, b);
					if (params->per_mask())
						kernels->halo_exchange(x);
				};
			}
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

		this->kreg->halo_setup(sop.grid(), &halo_ctx);
		sop.halo_ctx = halo_ctx;
		this->kreg->halo_stencil_exchange(sop);

		for (auto i :range<std::size_t>(this->levels.size()-1)) {
			this->levels.get(i+1).x.halo_ctx = halo_ctx;
			this->levels.get(i+1).b.halo_ctx = halo_ctx;
			this->levels.get(i+1).res.halo_ctx = halo_ctx;
			this->levels.get(i+1).A.halo_ctx = halo_ctx;
			this->levels.get(i+1).P.halo_ctx = halo_ctx;
		}
	}
	MPI_Comm comm;

private:
	bool cg_solver_lu;
	void *halo_ctx;
};

}}}

#endif
