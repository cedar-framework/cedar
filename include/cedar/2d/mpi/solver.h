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


	solver(mpi::stencil_op<fsten> & fop,
	       std::shared_ptr<config> conf,
	       service_manager<stypes> & outer_sman) : parent::multilevel(fop, conf), comm(fop.grid().comm)
	{
		this->kman = build_kernel_manager(*this->conf);

		// copy services (from outer redist solver) to new service manager
		service_manager<stypes> & newsman = kman->services();
		newsman.set_user_reg(outer_sman.get_user_reg());
		newsman.set<mpi::message_passing>(outer_sman.get_key<mpi::message_passing>());
		newsman.set<mpi::halo_exchange>(outer_sman.get_key<mpi::halo_exchange>());

		parent::setup(fop);
	}


	~solver() {}


	std::size_t compute_num_levels(mpi::stencil_op<fsten> & fop);


	virtual cdr2::mpi::grid_func solve(const cdr2::mpi::grid_func & b) override;

	virtual void solve(const cdr2::mpi::grid_func & b, cdr2::mpi::grid_func & x) override;

	virtual void setup_cg_solve() override;

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


	void setup_space(std::size_t nlevels);

	void setup_halo();

	void give_op(std::unique_ptr<stencil_op<fsten>> fop) {fop_ref = std::move(fop);}

	void apply_heirs(std::function<void(solver<nine_pt> &)> fun);

	MPI_Comm comm;


protected:
	std::unique_ptr<stencil_op<fsten>> fop_ref;
	std::shared_ptr<redist_solver<multilevel_wrapper<solver<nine_pt>>>> heir;
};

}}}

#endif
