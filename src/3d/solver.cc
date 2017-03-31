#include <algorithm>

// TODO: remove include
#include <cedar/3d/relax/setup_planes.h>
#include <cedar/3d/relax/relax_planes.h>

#include <cedar/3d/kernel/factory.h>
#include <cedar/3d/kernel/registry.h>
#include <cedar/3d/solver.h>

using namespace cedar;
using namespace cedar::cdr3;

void solver::setup_space(int nlevels)
{
	{
		grid_stencil & fsten = levels.back().A.stencil();
		levels.back().res = grid_func(fsten.shape(0), fsten.shape(1), fsten.shape(2));
	}

	for (auto i : range(nlevels - 1)) {
		auto & fop = levels[i].A;
		grid_stencil & sten = fop.stencil();
		auto nx = sten.shape(0);
		auto ny = sten.shape(1);
		auto nz = sten.shape(2);
		auto nxc = (nx-1)/2. + 1;
		auto nyc = (ny-1)/2. + 1;
		auto nzc = (nz-1)/2. + 1;

		auto cop = stencil_op(nxc, nyc, nzc);
		auto P = inter::prolong_op(nxc, nyc, nzc);
		std::array<relax_stencil, 2> SOR{{relax_stencil(nx, ny, nz),
					relax_stencil(nx, ny, nz)}};
		grid_stencil & st = cop.stencil();

		log::debug << "Created coarse grid with dimensions: " << st.shape(0)
		           << ", " << st.shape(1) << ", " << st.shape(2) << std::endl;

		levels.back().P = std::move(P);
		levels.back().SOR = std::move(SOR);
		cop.set_registry(kreg);

		levels.emplace_back(std::move(cop),inter::prolong_op());
		levels.back().x = grid_func(nxc, nyc, nzc);
		levels.back().res = grid_func(nxc, nyc, nzc);
		levels.back().b = grid_func(nxc, nyc, nzc);
	}

	auto & cop = levels.back().A;
	auto & cop_sten = cop.stencil();
	auto nxc = cop_sten.shape(0);
	auto nyc = cop_sten.shape(1);
	auto nzc = cop_sten.shape(2);
	ABD = grid_func(nxc*(nyc+1)+2, nxc*nyc*nzc, 0);
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
	float nxc, nyc, nzc;
	int ng = 0;
	int min_coarse = conf->get<int>("solver.min-coarse", 3);
	grid_stencil & sten = fop.stencil();

	auto nx = sten.shape(0);
	auto ny = sten.shape(1);
	auto nz = sten.shape(2);

	do {
		ng++;
		nxc = (nx-1)/(1<<ng) + 1;
		nyc = (ny-1)/(1<<ng) + 1;
		nzc = (nz-1)/(1<<ng) + 1;
	} while(std::min({nxc,nyc,nzc}) >= min_coarse);

	return ng;
}


void solver::setup_relax_plane(stencil_op & sop, bmg_level & level)
{
	kernel::impls::setup_relax_xy(sop, level.planes_xy);
	kernel::impls::setup_relax_xz(sop, level.planes_xz);
	kernel::impls::setup_relax_yz(sop, level.planes_yz);
}


void solver::relax_plane(const stencil_op & so, grid_func & x,
                         const grid_func & b, cycle::Dir cdir,
                         bmg_level & level)
{
	if (cdir == cycle::Dir::DOWN) {
		kernel::impls::relax_xy(so, x, b, cdir, level.planes_xy);
		kernel::impls::relax_yz(so, x, b, cdir, level.planes_yz);
		kernel::impls::relax_xz(so, x, b, cdir, level.planes_xz);
	} else {
		kernel::impls::relax_xz(so, x, b, cdir, level.planes_xz);
		kernel::impls::relax_yz(so, x, b, cdir, level.planes_yz);
		kernel::impls::relax_xy(so, x, b, cdir, level.planes_xy);
	}

}
