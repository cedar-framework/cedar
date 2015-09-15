#include <algorithm>

#include "boxmg-common.h"

#include "kernel/halo.h"
#include "kernel/mpi/factory.h"
#include "solver/mpi/boxmg.h"
#include "solver/boxmg.h"

using namespace boxmg;
using namespace boxmg::bmg2d;
using namespace boxmg::bmg2d::solver::mpi;

BoxMG::BoxMG(core::mpi::StencilOp&& fop) : comm(fop.grid().comm)
{
	kreg = kernel::mpi::factory::from_config(conf);
	fop.set_registry(kreg);

	levels.emplace_back(std::move(fop), inter::mpi::ProlongOp());

	auto num_levels = BoxMG::compute_num_levels(levels[0].A);
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
		levels[i].A.halo_ctx = halo_ctx;
		if (i != num_levels-1) {
			levels[i].P.halo_ctx = halo_ctx;
			std::array<core::RelaxStencil,2> SOR{{core::RelaxStencil(levels[i].A.stencil().shape(0), levels[i].A.stencil().shape(1)),core::RelaxStencil(levels[i].A.stencil().shape(0), levels[i].A.stencil().shape(1))}};
			levels[i].SOR = std::move(SOR);
			int kf = num_levels - i;
			kernels->setup_interp(kf,kf-1,num_levels,
			                     levels[i].A, levels[i+1].A,
			                     levels[i].P);
			kernels->galerkin_prod(kf, kf-1, num_levels, levels[i].P, levels[i].A, levels[i+1].A);
			auto relax_type = conf.get<std::string>("solver.relaxation", "point");

			if (relax_type == "point")
				kernels->setup_relax(levels[i].A,  levels[i].SOR[0]);
			else if (relax_type == "line-x")
				kernels->setup_relax_x(levels[i].A, levels[i].SOR[0]);
			else if (relax_type == "line-y")
				kernels->setup_relax_y(levels[i].A, levels[i].SOR[0]);
			else {
				kernels->setup_relax_x(levels[i].A, levels[i].SOR[0]);
				kernels->setup_relax_y(levels[i].A, levels[i].SOR[1]);
			}
			int nrelax_pre = conf.get<int>("solver.cycle.nrelax-pre", 2);
			int nrelax_post = conf.get<int>("solver.cycle.nrelax-post", 1);
			levels[i].presmoother = [&,i,nrelax_pre,kernels,relax_type](const core::DiscreteOp &A, core::GridFunc &x, const core::GridFunc &b) {
				const core::StencilOp & av = dynamic_cast<const core::StencilOp &>(A);
				for (auto j : range(nrelax_pre)) {
					(void)j;
					if (relax_type == "point")
						kernels->relax(av, x, b, levels[i].SOR[0], cycle::Dir::DOWN);
					else if (relax_type == "line-x")
						kernels->relax_lines_x(av, x, b, levels[i].SOR[0], cycle::Dir::DOWN);
					else if (relax_type == "line-y")
						kernels->relax_lines_y(av, x, b, levels[i].SOR[0], cycle::Dir::DOWN);
					else {
						kernels->relax_lines_x(av, x, b, levels[i].SOR[0], cycle::Dir::DOWN);
						kernels->relax_lines_y(av, x, b, levels[i].SOR[1], cycle::Dir::DOWN);
					}
				}
			};
			levels[i].postsmoother = [&,i,nrelax_post,kernels,relax_type](const core::DiscreteOp &A, core::GridFunc &x, const core::GridFunc&b) {

				const core::StencilOp & av = dynamic_cast<const core::StencilOp &>(A);
				for (auto j: range(nrelax_post)) {
					(void)j;
					if (relax_type == "point")
						kernels->relax(av, x, b, levels[i].SOR[0], cycle::Dir::UP);
					else if (relax_type == "line-x")
						kernels->relax_lines_x(av, x, b, levels[i].SOR[0], cycle::Dir::UP);
					else if (relax_type == "line-y")
						kernels->relax_lines_y(av, x, b, levels[i].SOR[0], cycle::Dir::UP);
					else {
						kernels->relax_lines_y(av, x, b, levels[i].SOR[1], cycle::Dir::UP);
						kernels->relax_lines_x(av, x, b, levels[i].SOR[0], cycle::Dir::UP);
					}
				}
			};
		}
	}

	auto & cop = levels.back().A;
	std::shared_ptr<solver::BoxMG> cg_bmg;
	kernels->setup_cg_boxmg(cop, &cg_bmg);

	coarse_solver = [&,cg_bmg,kernels](const core::DiscreteOp &A, core::mpi::GridFunc &x, const core::mpi::GridFunc &b) {
		const core::mpi::StencilOp &av = dynamic_cast<const core::mpi::StencilOp&>(A);
		auto &b_rw = const_cast<core::mpi::GridFunc&>(b);
		b_rw.halo_ctx = av.halo_ctx;
		kernels->solve_cg_boxmg(*cg_bmg, x, b);
		core::mpi::GridFunc residual = av.residual(x,b);
		log::info << "Level 0 residual norm: " << residual.lp_norm<2>() << std::endl;
	};
}


void BoxMG::add_level(core::mpi::StencilOp & fop, int num_levels)
{
	int kc = num_levels - levels.size() - 1;

	core::mpi::GridTopo & fgrid = fop.grid();
	auto cgrid = std::make_shared<core::mpi::GridTopo>(fgrid.get_igrd(), kc, num_levels);
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


	auto cop = core::mpi::StencilOp(cgrid);
	levels.back().P = inter::mpi::ProlongOp(cgrid);

	cop.set_registry(kreg);

	levels.emplace_back(std::move(cop),inter::mpi::ProlongOp());
}


int BoxMG::compute_num_levels(core::mpi::StencilOp & fop)
{
	int ng;
	auto min_coarse = conf.get<len_t>("solver.min_coarse", 3);

	auto kernels = kernel_registry();

	kernels->setup_nog(fop.grid(), min_coarse, &ng);

	return ng;
}


std::shared_ptr<kernel::Registry> BoxMG::kernel_registry()
{
	return std::static_pointer_cast<kernel::Registry>(kreg);
}


core::mpi::GridFunc BoxMG::solve(const core::mpi::GridFunc & b)
{
	auto kernels = kernel_registry();
	kernels->halo_exchange(b, halo_ctx);
	return MultiLevel<BoxMGLevel,core::mpi::GridFunc>::solve(b);
}
