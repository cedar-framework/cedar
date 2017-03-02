#include <algorithm>

#include <boxmg/types.h>

#include <boxmg/3d/mpi/halo.h>
#include <boxmg/3d/kernel/mpi/factory.h>
#include <boxmg/3d/mpi/solver.h>
#include <boxmg/perf/predict.h>
#include <boxmg/3d/solver.h>
#include <boxmg/3d/mpi/redist_solver.h>

using namespace boxmg;
using namespace boxmg::bmg3;


void mpi::solver::setup_space(int nlevels)
{
	levels.back().A.grid().grow(nlevels);
	levels[0].res = bmg3::mpi::grid_func(levels[0].A.grid_ptr());
	for (auto i : range(nlevels-1)) {
		auto & fop = levels.back().A;
		int kc = nlevels - levels.size() - 1;

		grid_topo & fgrid = fop.grid();
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


		auto cop = bmg3::mpi::stencil_op(cgrid);
		levels.back().P = inter::mpi::prolong_op(cgrid);
		std::array<bmg3::relax_stencil,2> SOR{{bmg3::relax_stencil(levels[i].A.stencil().shape(0), levels[i].A.stencil().shape(1), levels[i].A.stencil().shape(2)),bmg3::relax_stencil(levels[i].A.stencil().shape(0), levels[i].A.stencil().shape(1), levels[i].A.stencil().shape(2))}};
		levels[i].SOR = std::move(SOR);

		cop.set_registry(kreg);

		levels.emplace_back(std::move(cop),inter::mpi::prolong_op());
		levels.back().x = bmg3::mpi::grid_func(levels.back().A.grid_ptr());
		levels.back().b = bmg3::mpi::grid_func(levels.back().A.grid_ptr());
		levels.back().res = bmg3::mpi::grid_func(levels.back().A.grid_ptr());
	}

	setup_halo();
}


void mpi::solver::setup_cg_solve()
{
	auto & cop = levels.back().A;
	std::string cg_solver_str = conf->get<std::string>("solver.cg-solver", "LU");

	if (cg_solver_str == "LU" or cop.grid().nproc() == 1)
		cg_solver_lu = true;
	else
		cg_solver_lu = false;

	if (cg_solver_lu) {
		auto & coarse_topo = cop.grid();
		auto nxc = coarse_topo.nglobal(0);
		auto nyc = coarse_topo.nglobal(1);
		auto nzc = coarse_topo.nglobal(2);
		ABD = mpi::grid_func(nxc*(nyc+1)+2, nxc*nyc*nzc, 0);
		bbd = new real_t[ABD.len(1)];
		multilevel::setup_cg_solve();
	} else {
		auto kernels = kernel_registry();
		auto & fgrid = levels[0].A.grid();

		auto choice = choose_redist<3>(*conf,
		                               std::array<int, 3>({fgrid.nproc(0), fgrid.nproc(1), fgrid.nproc(2)}),
		                               std::array<len_t, 3>({fgrid.nglobal(0), fgrid.nglobal(1), fgrid.nglobal(2)}));

		MPI_Bcast(choice.data(), 3, MPI_INT, 0, fgrid.comm);
		log::status << "Redistributing to " << choice[0] << " x " << choice[1] << " x " << choice[2]
		            << " cores" << std::endl;

		std::shared_ptr<mpi::redist_solver> cg_bmg;
		auto cg_conf = conf->getconf("cg-config");
		if (!cg_conf)
			cg_conf = conf;
		std::vector<int> nblocks(choice.begin(), choice.end());
		kernels->setup_cg_redist(cop, cg_conf, &cg_bmg, nblocks);
		coarse_solver = [&,cg_bmg,kernels](const discrete_op<mpi::grid_func> &A, mpi::grid_func &x, const mpi::grid_func &b) {
			const bmg3::mpi::stencil_op &av = dynamic_cast<const bmg3::mpi::stencil_op&>(A);
			kernels->solve_cg_redist(*cg_bmg, x, b);
			// bmg3::mpi::grid_func residual = av.residual(x,b);
			// log::info << "Level 0 residual norm: " << residual.lp_norm<2>() << std::endl;
		};
	}
}


mpi::solver::solver(bmg3::mpi::stencil_op&& fop) : comm(fop.grid().comm)
{
	timer setup_timer("Setup");
	setup_timer.begin();

	kreg = kernel::mpi::factory::from_config(*conf);

	setup(std::move(fop));

	setup_timer.end();
}


mpi::solver::solver(bmg3::mpi::stencil_op&& fop,
                    std::shared_ptr<config::reader> cfg) : multilevel(cfg), comm(fop.grid().comm)
{
	timer setup_timer("Setup");
	setup_timer.begin();

	kreg = kernel::mpi::factory::from_config(*conf);

	setup(std::move(fop));

	setup_timer.end();
}


int mpi::solver::compute_num_levels(bmg3::mpi::stencil_op & fop)
{
	int ng;
	auto min_coarse = conf->get<len_t>("solver.min-coarse", 3);

	auto kernels = kernel_registry();

	kernels->setup_nog(fop.grid(), min_coarse, &ng);

	return ng;
}


bmg3::mpi::grid_func mpi::solver::solve(const bmg3::mpi::grid_func & b)
{
	auto kernels = kernel_registry();
	kernels->halo_exchange(b, halo_ctx);
	return multilevel::solve(b);
}


void mpi::solver::solve(const bmg3::mpi::grid_func & b, bmg3::mpi::grid_func & x)
{
	auto kernels = kernel_registry();
	kernels->halo_exchange(b, halo_ctx);
	return multilevel::solve(b, x);
}


void mpi::solver::setup_halo()
{
	auto & sop = levels[0].A;

	kreg->halo_setup(sop.grid(), &halo_ctx);
	sop.halo_ctx = halo_ctx;
	kreg->halo_stencil_exchange(sop);

	for (auto i :range(levels.size()-1)) {
		levels[i+1].x.halo_ctx = halo_ctx;
		levels[i+1].b.halo_ctx = halo_ctx;
		levels[i+1].res.halo_ctx = halo_ctx;
		levels[i+1].A.halo_ctx = halo_ctx;
		levels[i].P.halo_ctx = halo_ctx;
	}
}
