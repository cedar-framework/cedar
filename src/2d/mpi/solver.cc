#include <algorithm>

#include <boxmg/types.h>

#include <boxmg/2d/mpi/halo.h>
#include <boxmg/2d/kernel/mpi/factory.h>
#include <boxmg/2d/mpi/solver.h>
#include <boxmg/2d/solver.h>

using namespace boxmg;
using namespace boxmg::bmg2d;


void mpi::solver::setup_space(int nlevels)
{
	levels[0].res = bmg2d::mpi::grid_func(levels[0].A.grid_ptr());
	for (auto i : range(nlevels-1)) {
		auto & fop = levels.back().A;
		int kc = nlevels - levels.size() - 1;

		grid_topo & fgrid = fop.grid();
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


		auto cop = bmg2d::mpi::stencil_op(cgrid);
		levels.back().P = inter::mpi::prolong_op(cgrid);
		levels.back().P.halo_ctx = levels.back().A.halo_ctx;
		std::array<bmg2d::relax_stencil,2> SOR{{bmg2d::relax_stencil(levels[i].A.stencil().shape(0),
		                                                             levels[i].A.stencil().shape(1)),
					bmg2d::relax_stencil(levels[i].A.stencil().shape(0), levels[i].A.stencil().shape(1))}};
		levels.back().SOR = std::move(SOR);

		cop.set_registry(kreg);
		cop.halo_ctx = levels.back().A.halo_ctx;

		levels.emplace_back(std::move(cop),inter::mpi::prolong_op());
		levels.back().x = bmg2d::mpi::grid_func(levels.back().A.grid_ptr());
		levels.back().b = bmg2d::mpi::grid_func(levels.back().A.grid_ptr());
		levels.back().x.halo_ctx = levels.back().A.halo_ctx;
		levels.back().b.halo_ctx = levels.back().A.halo_ctx;
	}

	{
		std::string cg_solver_str = conf.get<std::string>("solver.cg-solver", "LU");
		if (cg_solver_str == "LU")
			cg_solver_lu = true;
		else
			cg_solver_lu = false;
	}


	if (cg_solver_lu) {
		auto & cop = levels.back().A;
		auto & coarse_topo = cop.grid();
		auto nxc = coarse_topo.nglobal(0);
		auto nyc = coarse_topo.nglobal(1);
		ABD = mpi::grid_func(nxc+2, nxc*nyc);
		bbd = new real_t[ABD.len(1)];
	}
}


void mpi::solver::setup_cg_solve()
{
	if (cg_solver_lu) {
		multilevel::setup_cg_solve();
	} else {
		auto kernels = kernel_registry();
		auto & cop = levels.back().A;
		std::shared_ptr<bmg2d::solver> cg_bmg;
		kernels->setup_cg_boxmg(cop, &cg_bmg);
		coarse_solver = [&,cg_bmg,kernels](const discrete_op<mpi::grid_func> &A, mpi::grid_func &x, const mpi::grid_func &b) {
			const bmg2d::mpi::stencil_op &av = dynamic_cast<const bmg2d::mpi::stencil_op&>(A);
			kernels->solve_cg_boxmg(*cg_bmg, x, b);
			bmg2d::mpi::grid_func residual = av.residual(x,b);
			log::info << "Level 0 residual norm: " << residual.lp_norm<2>() << std::endl;
		};
	}
}


mpi::solver::solver(bmg2d::mpi::stencil_op&& fop) : comm(fop.grid().comm)
{
	timer setup_timer("Setup");
	setup_timer.begin();

	kreg = kernel::mpi::factory::from_config(conf);

	kreg->halo_setup(levels[0].A.grid(), &halo_ctx);
	fop.halo_ctx = halo_ctx;
	kreg->halo_stencil_exchange(fop);

	setup(std::move(fop));

	setup_timer.end();
}


int mpi::solver::compute_num_levels(bmg2d::mpi::stencil_op & fop)
{
	int ng;
	auto min_coarse = conf.get<len_t>("solver.min-coarse", 3);

	auto kernels = kernel_registry();

	kernels->setup_nog(fop.grid(), min_coarse, &ng);

	return ng;
}


bmg2d::mpi::grid_func mpi::solver::solve(const bmg2d::mpi::grid_func & b)
{
	auto kernels = kernel_registry();
	kernels->halo_exchange(b, halo_ctx);
	return multilevel::solve(b);
}


void mpi::solver::solve(const bmg2d::mpi::grid_func & b, bmg2d::mpi::grid_func & x)
{
	auto kernels = kernel_registry();
	kernels->halo_exchange(b, halo_ctx);
	return multilevel::solve(b, x);
}