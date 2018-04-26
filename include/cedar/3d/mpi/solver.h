#ifndef CEDAR_3D_SOLVER_MPI_CEDAR_H
#define CEDAR_3D_SOLVER_MPI_CEDAR_H

#include <algorithm>
#include <memory>
#include <mpi.h>
#include <array>

#include <cedar/multilevel.h>
#include <cedar/level.h>
#include <cedar/perf/predict.h>
#include <cedar/3d/level_container.h>
#include <cedar/3d/kernel_manager.h>
#include <cedar/3d/mpi/types.h>
#include <cedar/3d/redist/multilevel_wrapper.h>
#include <cedar/3d/redist/cholesky_solver.h>
#include <cedar/3d/redist/setup_redist.h>

namespace cedar { namespace cdr3 { namespace mpi {


template<class sten>
	struct level3mpi : public level<sten, stypes>
{
	using parent = level<sten, stypes>;
level3mpi(topo_ptr topo ): parent::level(topo)
	{
		this->SOR = {{relax_stencil(topo->nlocal(0) - 2, topo->nlocal(1) - 2, topo->nlocal(2) - 2),
		              relax_stencil(topo->nlocal(0) - 2, topo->nlocal(1) - 2, topo->nlocal(2) - 2)}};
		this->R.associate(&this->P);
	}

level3mpi(stencil_op<sten> & A) : parent::level(A)
	{
		this->res = mpi::grid_func(A.grid_ptr());
		this->SOR = {{relax_stencil(A.shape(0), A.shape(1), A.shape(2)),
		              relax_stencil(A.shape(0), A.shape(1), A.shape(2)),}};
	}
};

template<class fsten>
class solver: public multilevel<exec_mode::mpi, level_container<level3mpi, fsten>,
                                fsten, cdr3::mpi::solver<fsten>>
{
public:
	using parent = multilevel<exec_mode::mpi, level_container<level3mpi, fsten>,
	                          fsten, cdr3::mpi::solver<fsten>>;
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
	~solver() {}

	std::size_t compute_num_levels(cdr3::mpi::stencil_op<fsten> & fop)
	{
		int ng;
		auto min_coarse = this->settings.min_coarse;

		kman->template run<setup_nog>(fop.grid(), min_coarse, &ng);

		return ng;
	}


	virtual cdr3::mpi::grid_func solve(const cdr3::mpi::grid_func &b) override
	{
		auto & bd = const_cast<grid_func&>(b);
		kman->template run<halo_exchange>(bd);
		return parent::solve(b);
	}


	virtual void solve(const cdr3::mpi::grid_func &b, cdr3::mpi::grid_func &x) override
	{
		auto & bd = const_cast<grid_func&>(b);
		kman->template run<halo_exchange>(bd);
		return parent::solve(b, x);
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
		for (auto i : range<std::size_t>(nlevels - 1)) {
			auto & fgrid = this->get_grid(i);
			if (i == 0)
				fgrid.grow(nlevels);

			int kc = nlevels - i - 2;

			auto cgrid = std::make_shared<grid_topo>(fgrid.get_igrd(), kc, nlevels);
			cgrid->comm = fgrid.comm;

			len_t NLxg = fgrid.nlocal(0) - 2;
			len_t NLyg = fgrid.nlocal(1) - 2;
			len_t NLzg = fgrid.nlocal(2) - 2;
			len_t NGxg = (fgrid.nglobal(0) - 1) / 2 + 2;
			len_t NGyg = (fgrid.nglobal(1) - 1) / 2 + 2;
			len_t NGzg = (fgrid.nglobal(2) - 1) / 2 + 2;

			cgrid->nglobal(0) = NGxg;
			cgrid->nglobal(1) = NGyg;
			cgrid->nglobal(2) = NGzg;

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

			if (fgrid.is(2) % 2 == 1) {
				cgrid->is(2) = (fgrid.is(2)+1) / 2;
				NLzg = (NLzg+1) / 2;
			} else {
				cgrid->is(2) = fgrid.is(2) / 2 + 1;
				if (NLzg % 2 == 1) NLzg = (NLzg - 1) / 2;
				else NLzg = (NLzg+1)/2;
			}

			cgrid->nlocal(0) = NLxg + 2;
			cgrid->nlocal(1) = NLyg + 2;
			cgrid->nlocal(2) = NLzg + 2;

			cgrid->nproc(0) = fgrid.nproc(0);
			cgrid->nproc(1) = fgrid.nproc(1);
			cgrid->nproc(2) = fgrid.nproc(2);
			cgrid->coord(0) = fgrid.coord(0);
			cgrid->coord(1) = fgrid.coord(1);
			cgrid->coord(2) = fgrid.coord(2);

			this->levels.add(cgrid);
		}

		setup_halo();
	}


	virtual void setup_cg_solve() override
	{
		auto params = build_kernel_params(*(this->conf));
		auto & cop = this->levels.get(this->levels.size() - 1).A;

		auto cg_conf = this->settings.coarse_config;

		if (this->settings.coarse_solver == ml_settings::cg_type::redist) {
			auto & fgrid = cop.grid();

			auto choice = choose_redist<3>(*this->conf,
			                               std::array<int, 3>({fgrid.nproc(0), fgrid.nproc(1), fgrid.nproc(2)}),
			                               std::array<len_t, 3>({fgrid.nglobal(0), fgrid.nglobal(1), fgrid.nglobal(2)}));

			MPI_Bcast(choice.data(), 3, MPI_INT, 0, fgrid.comm);
			if ((choice[0] != 1) or (choice[1] != 1) or (choice[2] != 1)) {
				log::status << "Redistributing to " << choice[0] << " x " << choice[1] << " x " << choice[2]
				            << " cores" << std::endl;
				using inner_solver = multilevel_wrapper<mpi::solver<xxvii_pt>>;
				this->coarse_solver = create_redist_solver<inner_solver>(kman,
				                                                         *this->conf,
				                                                         cop,
				                                                         cg_conf,
				                                                         choice);
				return;
			}
		}

		std::array<int, 3> choice {{1,1,1}};
		log::status << "Redistributing to " << choice[0] << " x " << choice[1] << " x " << choice[2]
		            << " cores" << std::endl;
		if (this->settings.coarse_solver == ml_settings::cg_type::lu) {
			this->coarse_solver = create_redist_solver<cholesky_solver>(kman,
			                                                            *this->conf,
			                                                            cop,
			                                                            cg_conf,
			                                                            choice);
		} else {
			using inner_solver = multilevel_wrapper<cdr3::solver<xxvii_pt>>;
			this->coarse_solver = create_redist_solver<inner_solver>(kman,
			                                                         *this->conf,
			                                                         cop,
			                                                         cg_conf,
			                                                         choice);
		}
	}


	void setup_halo()
	{
		auto & sop = this->levels.template get<fsten>(0).A;

		std::vector<topo_ptr> topos;
		topos.push_back(sop.grid_ptr());

		for (std::size_t i = 1; i < this->nlevels(); i++)
			topos.push_back(this->levels.get(i).A.grid_ptr());

		this->kman->template setup<halo_exchange>(topos);
		this->kman->template run<halo_exchange>(sop);
	}

	MPI_Comm comm;

private:
	void *halo_ctx;
};

}}}

#endif
