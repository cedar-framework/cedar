#include <algorithm>

#include <cedar/types.h>

#include <cedar/2d/mpi/halo.h>
#include <cedar/2d/kernel/mpi/factory.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/perf/predict.h>
#include <cedar/2d/solver.h>

using namespace cedar;
using namespace cedar::cdr2;


void mpi::solver::setup_space(int nlevels)
{
	levels.back().A.grid().grow(nlevels);
	levels[0].res = cdr2::mpi::grid_func(levels[0].A.grid_ptr());
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


		auto cop = cdr2::mpi::stencil_op(cgrid);
		levels.back().P = inter::mpi::prolong_op(cgrid);
		std::array<cdr2::relax_stencil,2> SOR{{cdr2::relax_stencil(levels[i].A.stencil().shape(0),
		                                                             levels[i].A.stencil().shape(1)),
					cdr2::relax_stencil(levels[i].A.stencil().shape(0), levels[i].A.stencil().shape(1))}};
		levels.back().SOR = std::move(SOR);

		cop.set_registry(kreg);

		levels.emplace_back(std::move(cop),inter::mpi::prolong_op());
		levels.back().x = cdr2::mpi::grid_func(levels.back().A.grid_ptr());
		levels.back().b = cdr2::mpi::grid_func(levels.back().A.grid_ptr());
		levels.back().res = cdr2::mpi::grid_func(levels.back().A.grid_ptr());
	}

	setup_halo();
}


void mpi::solver::setup_cg_solve()
{
	auto & cop = levels.back().A;
	std::string cg_solver_str = conf->get<std::string>("solver.cg-solver", "LU");

	if (cg_solver_str == "LU")
		cg_solver_lu = true;
	else
		cg_solver_lu = false;

	if (cg_solver_lu) {
		auto & cop = levels.back().A;
		auto & coarse_topo = cop.grid();
		auto nxc = coarse_topo.nglobal(0);
		auto nyc = coarse_topo.nglobal(1);
		ABD = mpi::grid_func(nxc+2, nxc*nyc);
		bbd = new real_t[ABD.len(1)];
		multilevel::setup_cg_solve();
	} else {
		auto kernels = kernel_registry();
		std::vector<int> nblocks({1,1});
		if (cg_solver_str == "redist") {
			auto & fgrid = levels[0].A.grid();
			{
				int rank;
				MPI_Comm_rank(fgrid.comm, &rank);
				nblocks = std::move(
					predict_redist(*conf, fgrid.nproc(0), fgrid.nproc(1), fgrid.nglobal(0), fgrid.nglobal(1))
					);
				MPI_Bcast(nblocks.data(), 2, MPI_INT, 0, fgrid.comm);
				log::status << "Redistributing to " << nblocks[0] << " x " << nblocks[1] << " cores" << std::endl;
			}
		}
		if (cg_solver_str == "boxmg" or nblocks[0]*nblocks[1] == 1) {
			std::shared_ptr<cdr2::solver> cg_bmg;
			auto cg_conf = conf->getconf("cg-config");
			if (!cg_conf)
				cg_conf = conf;
			kernels->setup_cg_boxmg(cop, cg_conf, &cg_bmg);
			coarse_solver = [&,cg_bmg,kernels](const discrete_op<mpi::grid_func> &A, mpi::grid_func &x, const mpi::grid_func &b) {
				const cdr2::mpi::stencil_op &av = dynamic_cast<const cdr2::mpi::stencil_op&>(A);
				kernels->solve_cg_boxmg(*cg_bmg, x, b);
				cdr2::mpi::grid_func residual = av.residual(x,b);
				log::info << "Level 0 residual norm: " << residual.lp_norm<2>() << std::endl;
			};
		} else if (cg_solver_str == "redist") {
			std::shared_ptr<mpi::redist_solver> cg_bmg;
			auto cg_conf = conf->getconf("cg-config");
			if (!cg_conf)
				cg_conf = conf;
			kernels->setup_cg_redist(cop, cg_conf, &cg_bmg, nblocks);
			coarse_solver = [&,cg_bmg,kernels](const discrete_op<mpi::grid_func> &A, mpi::grid_func &x, const mpi::grid_func &b) {
				const cdr2::mpi::stencil_op &av = dynamic_cast<const cdr2::mpi::stencil_op&>(A);
				kernels->solve_cg_redist(*cg_bmg, x, b);
				cdr2::mpi::grid_func residual = av.residual(x,b);
				log::info << "Level 0 residual norm: " << residual.lp_norm<2>() << std::endl;
			};
		}
	}
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




mpi::solver::solver(cdr2::mpi::stencil_op&& fop) : comm(fop.grid().comm)
{
	kreg = kernel::mpi::factory::from_config(*conf);

	setup(std::move(fop));
}


mpi::solver::solver(cdr2::mpi::stencil_op&& fop,
                    std::shared_ptr<config::reader> cfg) : multilevel(cfg), comm(fop.grid().comm)
{
	kreg = kernel::mpi::factory::from_config(*conf);

	setup(std::move(fop));
}


int mpi::solver::compute_num_levels(cdr2::mpi::stencil_op & fop)
{
	int ng;
	auto min_coarse = conf->get<len_t>("solver.min-coarse", 3);

	auto kernels = kernel_registry();

	kernels->setup_nog(fop.grid(), min_coarse, &ng);

	return ng;
}


cdr2::mpi::grid_func mpi::solver::solve(const cdr2::mpi::grid_func & b)
{
	auto kernels = kernel_registry();
	kernels->halo_exchange(b, halo_ctx);
	return multilevel::solve(b);
}


void mpi::solver::solve(const cdr2::mpi::grid_func & b, cdr2::mpi::grid_func & x)
{
	auto kernels = kernel_registry();
	kernels->halo_exchange(b, halo_ctx);
	return multilevel::solve(b, x);
}
