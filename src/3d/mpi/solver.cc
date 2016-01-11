#include <algorithm>

#include <boxmg/types.h>

#include <boxmg/3d/mpi/halo.h>
#include <boxmg/3d/kernel/mpi/factory.h>
#include <boxmg/3d/mpi/solver.h>
#include <boxmg/3d/solver.h>

using namespace boxmg;
using namespace boxmg::bmg3;

mpi::solver::solver(bmg3::mpi::stencil_op&& fop) : comm(fop.grid().comm)
{
	Timer setup_timer("Setup");
	setup_timer.begin();

	kreg = kernel::mpi::factory::from_config(conf);
	fop.set_registry(kreg);

	levels.emplace_back(std::move(fop), inter::mpi::prolong_op());

	auto num_levels = mpi::solver::compute_num_levels(levels[0].A);
	log::debug << "Using a " << num_levels << " level heirarchy" << std::endl;
	levels.reserve(num_levels);
	levels.back().A.grid().grow(num_levels);
	for (auto i: range(num_levels-1)) {
		add_level(levels.back().A, num_levels);
		levels[i].R.associate(&levels[i].P);
		levels[i].P.fine_op = &levels[i].A;
		levels[i].R.set_registry(kreg);
		levels[i].P.set_registry(kreg);
	}

	auto kernels = kernel_registry();
	kernels->halo_setup(levels[0].A.grid(), &halo_ctx);
	levels[0].A.halo_ctx = halo_ctx;
	kernels->halo_stencil_exchange(levels[0].A);

	for (auto i: range(num_levels)) {
		levels[i].res = bmg3::mpi::grid_func(levels[i].A.grid_ptr());
		if (i) {
			levels[i].x = bmg3::mpi::grid_func(levels[i].A.grid_ptr());
			levels[i].b = bmg3::mpi::grid_func(levels[i].A.grid_ptr());
		}
		levels[i].A.halo_ctx = halo_ctx;
		if (i != num_levels-1) {
			levels[i].P.halo_ctx = halo_ctx;
			std::array<bmg3::relax_stencil,2> SOR{{bmg3::relax_stencil(levels[i].A.stencil().shape(0), levels[i].A.stencil().shape(1), levels[i].A.stencil().shape(2)),bmg3::relax_stencil(levels[i].A.stencil().shape(0), levels[i].A.stencil().shape(1), levels[i].A.stencil().shape(2))}};
			levels[i].SOR = std::move(SOR);
			int kf = num_levels - i;
			kernels->setup_interp(kf,kf-1,num_levels,
			                     levels[i].A, levels[i+1].A,
			                     levels[i].P);
			kernels->galerkin_prod(kf, kf-1, num_levels, levels[i].P, levels[i].A, levels[i+1].A);
			auto relax_type = conf.get<std::string>("solver.relaxation", "point");

			if (relax_type == "point")
				kernels->setup_relax(levels[i].A,  levels[i].SOR[0]);
			else {
				log::error << "Plane relaxation not yet added" << std::endl;
			}
			int nrelax_pre = conf.get<int>("solver.cycle.nrelax-pre", 2);
			int nrelax_post = conf.get<int>("solver.cycle.nrelax-post", 1);
			levels[i].presmoother = [&,i,nrelax_pre,kernels,relax_type](const discrete_op<mpi::grid_func> &A, mpi::grid_func &x, const mpi::grid_func &b) {
				const bmg3::mpi::stencil_op & av = dynamic_cast<const bmg3::mpi::stencil_op &>(A);
				for (auto j : range(nrelax_pre)) {
					(void)j;
					if (relax_type == "point")
						kernels->relax(av, x, b, levels[i].SOR[0], cycle::Dir::DOWN);
					else {
						log::error << "Plane relaxation not yet added" << std::endl;
					}
				}
			};
			levels[i].postsmoother = [&,i,nrelax_post,kernels,relax_type](const discrete_op<mpi::grid_func> &A, mpi::grid_func &x, const mpi::grid_func&b) {

				const bmg3::mpi::stencil_op & av = dynamic_cast<const bmg3::mpi::stencil_op &>(A);
				for (auto j: range(nrelax_post)) {
					(void)j;
					if (relax_type == "point")
						kernels->relax(av, x, b, levels[i].SOR[0], cycle::Dir::UP);
					else {
						log::error << "Plane relaxation not yet added" << std::endl;
					}
				}
			};
		}
	}

	auto & cop = levels.back().A;
	{
		std::string cg_solver_str = conf.get<std::string>("solver.cg-solver", "LU");
		if (cg_solver_str == "LU")
			cg_solver_lu = true;
		else
			cg_solver_lu = false;
	}

	std::shared_ptr<bmg3::solver> cg_bmg;
	if (cg_solver_lu) {
		auto & coarse_topo = cop.grid();
		auto nxc = coarse_topo.nglobal(0);
		auto nyc = coarse_topo.nglobal(1);
		auto nzc = coarse_topo.nglobal(2);
		ABD = mpi::grid_func(nxc*(nyc+1)+2, nxc*nyc*nzc, 0);
		bbd = new real_t[ABD.len(1)];
		kernels->setup_cg_lu(cop, ABD);
	} else {
		kernels->setup_cg_boxmg(cop, &cg_bmg);
	}

	coarse_solver = [&,cg_bmg,kernels](const discrete_op<mpi::grid_func> &A, mpi::grid_func &x, const mpi::grid_func &b) {
		const bmg3::mpi::stencil_op &av = dynamic_cast<const bmg3::mpi::stencil_op&>(A);
		auto &b_rw = const_cast<bmg3::mpi::grid_func&>(b);
		b_rw.halo_ctx = av.halo_ctx;
		if (cg_solver_lu)
			kernels->solve_cg(x, b, ABD, bbd);
		else
			kernels->solve_cg_boxmg(*cg_bmg, x, b);
		bmg3::mpi::grid_func residual = av.residual(x,b);
		log::info << "Level 0 residual norm: " << residual.lp_norm<2>() << std::endl;
	};

	setup_timer.end();
}


void mpi::solver::add_level(bmg3::mpi::stencil_op & fop, int num_levels)
{
	int kc = num_levels - levels.size() - 1;

	grid_topo & fgrid = fop.grid();
	auto cgrid = std::make_shared<grid_topo>(fgrid.get_igrd(), kc, num_levels);
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

	cop.set_registry(kreg);

	levels.emplace_back(std::move(cop),inter::mpi::prolong_op());
}


int mpi::solver::compute_num_levels(bmg3::mpi::stencil_op & fop)
{
	int ng;
	auto min_coarse = conf.get<len_t>("solver.min-coarse", 3);

	auto kernels = kernel_registry();

	kernels->setup_nog(fop.grid(), min_coarse, &ng);

	return ng;
}


std::shared_ptr<boxmg::bmg3::kernel::mpi::registry> mpi::solver::kernel_registry()
{
	return std::static_pointer_cast<boxmg::bmg3::kernel::mpi::registry>(kreg);
}


bmg3::mpi::grid_func mpi::solver::solve(const bmg3::mpi::grid_func & b)
{
	auto kernels = kernel_registry();
	kernels->halo_exchange(b, halo_ctx);
	return multilevel<BoxMGLevel,bmg3::mpi::grid_func, boxmg::bmg3::kernel::mpi::registry>::solve(b);
}


void mpi::solver::solve(const bmg3::mpi::grid_func & b, bmg3::mpi::grid_func & x)
{
	auto kernels = kernel_registry();
	kernels->halo_exchange(b, halo_ctx);
	return multilevel<BoxMGLevel,bmg3::mpi::grid_func,boxmg::bmg3::kernel::mpi::registry>::solve(b, x);
}
