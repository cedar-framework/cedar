#include <algorithm>

#include "boxmg-common.h"

#include "kernel/calls.h"
#include "kernel/factory.h"
#include "solver/boxmg.h"

using namespace boxmg;
using namespace boxmg::bmg2d;
using namespace boxmg::bmg2d::solver;

BoxMG::BoxMG(core::StencilOp&& fop)
{
	kreg = kernel::factory::from_config(conf);
	levels.emplace_back(std::move(fop), inter::ProlongOp());
	levels.back().A.set_registry(kreg);
	{
		core::GridStencil & fsten = levels.back().A.stencil();
		levels.back().res = core::GridFunc(fsten.shape(0), fsten.shape(1));
	}

	auto num_levels = BoxMG::compute_num_levels(levels[0].A);
	log::debug << "Using a " << num_levels << " level heirarchy" << std::endl;
	levels.reserve(num_levels);
	for (auto i: range(num_levels-1)) {
		add_level(levels.back().A, num_levels);
		levels[i].R.associate(&levels[i].P);
		levels[i].P.fine_op = &levels[i].A;
		levels[i].R.set_registry(kreg);
		levels[i].P.set_registry(kreg);
	}

	auto kernels = kernel_registry();

	auto & cop = levels.back().A;
	auto & cop_sten = cop.stencil();
	auto nxc = cop_sten.shape(0);
	auto nyc = cop_sten.shape(1);
	ABD = core::GridFunc(nxc+2, nxc*nyc, 0);
	bbd = new real_t[ABD.len(1)];
	kernels->setup_cg_lu(cop, ABD);
	coarse_solver = [&,kernels](const core::DiscreteOp &A, core::GridFunc &x, const core::GridFunc &b) {
		kernels->solve_cg(x, b, ABD, bbd);
		const core::StencilOp &av = dynamic_cast<const core::StencilOp&>(A);
		core::GridFunc & residual = levels[levels.size()-1].res;
		av.residual(x,b,residual);
		log::info << "Level 0 residual norm: " << residual.lp_norm<2>() << std::endl;
	};
}


void BoxMG::add_level(core::StencilOp & fop, int num_levels)
{
	core::GridStencil & sten = fop.stencil();
	auto kernels = kernel_registry();
	auto nx = sten.shape(0);
	auto ny = sten.shape(1);
	auto nxc = (nx-1)/2. + 1;
	auto nyc = (ny-1)/2. + 1;
	int kc = num_levels - levels.size() - 1;
	int kf = kc + 1;

	auto cop = core::StencilOp(nxc, nyc);
	auto P = inter::ProlongOp(nxc, nyc);
	std::array<core::RelaxStencil, 2> SOR{{core::RelaxStencil(nx, ny),
				core::RelaxStencil(nx, ny)}};
	core::GridStencil & st = cop.stencil();

	log::debug << "Created coarse grid with dimensions: " << st.shape(0)
	          << ", " << st.shape(1) << std::endl;

	std::string relax_type = conf.get<std::string>("solver.relaxation", "point");

	kernels->setup_interp(kf, kc, num_levels, fop, cop, P);
	kernels->galerkin_prod(kf, kc, num_levels, P, fop, cop);

	if (relax_type == "point")
		kernels->setup_relax(fop, SOR[0]);
	else if (relax_type == "line-x")
		kernels->setup_relax_x(fop, SOR[0]);
	else if (relax_type == "line-y")
		kernels->setup_relax_y(fop, SOR[0]);
	else { // line-xy
		kernels->setup_relax_x(fop, SOR[0]);
		kernels->setup_relax_y(fop, SOR[1]);
	}

	levels.back().P = std::move(P);
	levels.back().SOR = std::move(SOR);

	auto lvl = levels.size() - 1;
	int nrelax_pre = conf.get<int>("solver.cycle.nrelax-pre", 2);
	int nrelax_post = conf.get<int>("solver.cycle.nrelax-post", 1);

	levels.back().presmoother = [&,lvl,nrelax_pre,kernels,relax_type](const core::DiscreteOp &A, core::GridFunc &x, const core::GridFunc&b) {
		const core::StencilOp & av = dynamic_cast<const core::StencilOp &>(A);
		for (auto i : range(nrelax_pre)) {
			(void) i;
			if (relax_type == "point")
				kernels->relax(av, x, b, levels[lvl].SOR[0], cycle::Dir::DOWN);
			else if (relax_type == "line-x")
				kernels->relax_lines_x(av, x, b, levels[lvl].SOR[0], cycle::Dir::DOWN);
			else if (relax_type == "line-y")
				kernels->relax_lines_y(av, x, b, levels[lvl].SOR[0], cycle::Dir::DOWN);
			else {
				kernels->relax_lines_x(av, x, b, levels[lvl].SOR[0], cycle::Dir::DOWN);
				kernels->relax_lines_y(av, x, b, levels[lvl].SOR[1], cycle::Dir::DOWN);
			}
		}
	};
	levels.back().postsmoother = [&,lvl,nrelax_post,kernels,relax_type](const core::DiscreteOp &A, core::GridFunc &x, const core::GridFunc&b) {

		const core::StencilOp & av = dynamic_cast<const core::StencilOp &>(A);
		for (auto i: range(nrelax_post)) {
			(void) i;
			if (relax_type == "point")
				kernels->relax(av, x, b, levels[lvl].SOR[0], cycle::Dir::DOWN);
			else if (relax_type == "line-x")
				kernels->relax_lines_x(av, x, b, levels[lvl].SOR[0], cycle::Dir::DOWN);
			else if (relax_type == "line-y")
				kernels->relax_lines_y(av, x, b, levels[lvl].SOR[0], cycle::Dir::DOWN);
			else {
				kernels->relax_lines_y(av, x, b, levels[lvl].SOR[1], cycle::Dir::DOWN);
				kernels->relax_lines_x(av, x, b, levels[lvl].SOR[0], cycle::Dir::DOWN);
			}
		}
	};

	cop.set_registry(kreg);

	levels.emplace_back(std::move(cop),inter::ProlongOp());
	levels.back().x = core::GridFunc(nxc, nyc);
	levels.back().res = core::GridFunc(nxc, nyc);
	levels.back().b = core::GridFunc(nxc, nyc);
}


int BoxMG::compute_num_levels(core::StencilOp & fop)
{
	float nxc, nyc;
	int ng = 0;
	int min_coarse = conf.get<int>("solver.min_coarse", 3);
	core::GridStencil & sten = fop.stencil();

	auto nx = sten.shape(0);
	auto ny = sten.shape(1);

	do {
		ng++;
		nxc = (nx-1)/(1<<ng) + 1;
		nyc = (ny-1)/(1<<ng) + 1;
	} while(std::min(nxc,nyc) >= min_coarse);

	return ng;
}


std::shared_ptr<kernel::Registry> BoxMG::kernel_registry()
{
	return std::static_pointer_cast<kernel::Registry>(kreg);
}
