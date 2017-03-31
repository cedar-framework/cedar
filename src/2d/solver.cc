#include <algorithm>

#include <cedar/2d/kernel/factory.h>
#include <cedar/2d/solver.h>
#include <cedar/2d/kernel/registry.h>

using namespace cedar;
using namespace cedar::cdr2;


void solver::setup_space(int nlevels)
{
	{
		grid_stencil & fsten = levels.back().A.stencil();
		levels.back().res = grid_func(fsten.shape(0), fsten.shape(1));
	}

	for (auto i : range(nlevels-1)) {
		auto & fop = levels[i].A;
		auto & sten = fop.stencil();
		auto nx = sten.shape(0);
		auto ny = sten.shape(1);
		auto nxc = (nx-1)/2. + 1;
		auto nyc = (ny-1)/2. + 1;

		auto cop = stencil_op(nxc, nyc);
		auto P = inter::prolong_op(nxc, nyc);
		std::array<relax_stencil, 2> SOR{{relax_stencil(nx, ny),
					relax_stencil(nx, ny)}};
		grid_stencil & st = cop.stencil();

		log::debug << "Created coarse grid with dimensions: " << st.shape(0)
		           << ", " << st.shape(1) << std::endl;

		levels.back().P = std::move(P);
		levels.back().SOR = std::move(SOR);

		cop.set_registry(kreg);
		levels.emplace_back(std::move(cop),inter::prolong_op());
		levels.back().x = grid_func(nxc, nyc);
		levels.back().res = grid_func(nxc, nyc);
		levels.back().b = grid_func(nxc, nyc);
	}

	auto & cop = levels.back().A;
	auto & cop_sten = cop.stencil();
	auto nxc = cop_sten.shape(0);
	auto nyc = cop_sten.shape(1);
	ABD = grid_func(nxc+2, nxc*nyc, 0);
	bbd = new real_t[ABD.len(1)];
}


solver::solver(stencil_op&& fop)
{
	kreg = kernel::factory::from_config(*conf);

	setup(std::move(fop));
}


solver::solver(stencil_op&& fop,
               std::shared_ptr<config::reader> cfg): multilevel(cfg)
{
	kreg = kernel::factory::from_config(*conf);

	setup(std::move(fop));
}


solver::~solver()
{
	delete[] bbd;
}


int solver::compute_num_levels(stencil_op & fop)
{
	float nxc, nyc;
	int ng = 0;
	int min_coarse = conf->get<int>("solver.min-coarse", 3);
	grid_stencil & sten = fop.stencil();

	auto nx = sten.shape(0);
	auto ny = sten.shape(1);

	do {
		ng++;
		nxc = (nx-1)/(1<<ng) + 1;
		nyc = (ny-1)/(1<<ng) + 1;
	} while(std::min(nxc,nyc) >= min_coarse);

	return ng;
}
