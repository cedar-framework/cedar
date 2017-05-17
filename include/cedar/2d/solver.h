#ifndef CEDAR_2D_CORE_SOLVER_H
#define CEDAR_2D_CORE_SOLVER_H

#include <array>
#include <functional>
#include <algorithm>

#include <cedar/multilevel.h>
#include <cedar/2d/stencil_op.h>
#include <cedar/2d/relax_stencil.h>
#include <cedar/2d/inter/prolong_op.h>
#include <cedar/2d/inter/restrict_op.h>
#include <cedar/2d/kernel/registry.h>

namespace cedar { namespace cdr2 {

template<class sten>
struct level
{
level(len_t nx, len_t ny) : Adata(nx, ny), A(Adata),
		P(nx, ny), x(nx, ny), res(nx, ny), b(nx, ny),
		SOR{{relax_stencil(A.shape(0), A.shape(1)),
				relax_stencil(A.shape(0), A.shape(1))}} {
	R.associate(&P);
}
level(stencil_op<sten> & A) : A(A),
		res(A.shape(0), A.shape(1)),
		SOR{{relax_stencil(A.shape(0), A.shape(1)),
				relax_stencil(A.shape(0), A.shape(1))}} {}
	stencil_op<sten> Adata;
	stencil_op<sten> & A;
	inter::prolong_op  P;
	inter::restrict_op R;
	grid_func x;
	grid_func res;
	grid_func b;
	std::array<relax_stencil, 2> SOR;

	std::function<void(const stencil_op<sten> & A, grid_func & x, const grid_func & b)> presmoother;
	std::function<void(const stencil_op<sten> & A, grid_func & x, const grid_func & b)> postsmoother;
};

template<class sten>
class level_container
{
public:
	template<class rsten> using value_type = level<rsten>;
	level_container(stencil_op<sten> & fine_op);
	void add(len_t nx, len_t ny) {
		lvls_nine.emplace_back(nx, ny);
	}
	template<class rsten=nine_pt> level<rsten>&  get(std::size_t i);
	std::size_t size() { return lvls_nine.size() + lvls_five.size(); }

protected:
	std::vector<level<nine_pt>> lvls_nine;
	std::vector<level<five_pt>> lvls_five;
};

template<> template<> inline level<nine_pt>& level_container<five_pt>::get<nine_pt>(std::size_t i)
{
	if (i==0) log::error << "fine grid operator is five point (not nine)!" << std::endl;
	#ifdef BOUNDS_CHECK
	return lvls_nine.at(i - lvls_five.size());
	#else
	return lvls_nine[i - lvls_five.size()];
	#endif
}

template<> template<> inline level<five_pt>& level_container<five_pt>::get<five_pt>(std::size_t i)
{
	if (i != 0) log::error << "coarse operators are nine point (not five)!" << std::endl;
	#ifdef BOUNDS_CHECK
	return lvls_five.at(0);
	#else
	return lvls_five[0];
	#endif
}

template<> template<> inline level<nine_pt>& level_container<nine_pt>::get<nine_pt>(std::size_t i)
{
	#ifdef BOUNDS_CHECK
	return lvls_nine.at(i);
	#else
	return lvls_nine[i];
	#endif
}

template<class fsten>
class solver: public multilevel<level_container, stencil_op, grid_func, kernel::registry, fsten, solver<fsten>>
{
public:
solver(stencil_op<fsten> & fop) : multilevel<level_container, stencil_op, grid_func, kernel::registry, fsten, solver<fsten>>(fop)
	{
		kreg = std::make_shared<kernel::registry>(*conf);
		multilevel<level_container, stencil_op, grid_func, kernel::registry, fsten, solver<fsten>>::setup();
	}
	solver(stencil_op<fsten> & fop,
	       std::shared_ptr<config::reader> conf) :
	multilevel<level_container, stencil_op, grid_func, kernel::registry, fsten, solver<fsten>>(fop, conf)
	{
		kreg = std::make_shared<kernel::registry>(*conf);
		multilevel<level_container, stencil_op, grid_func, kernel::registry, fsten, solver<fsten>>::setup();
	}
	~solver() { delete[] this->bbd; }
	int compute_num_levels(stencil_op<sten> & fop)
	{
		float nxc, nyc;
		int ng = 0;
		int min_coarse = conf->get<int>("solver.min-coarse", 3);

		auto nx = fop.shape(0);
		auto ny = fop.shape(1);

		do {
			ng++;
			nxc = (nx-1)/(1<<ng) + 1;
			nyc = (ny-1)/(1<<ng) + 1;
		} while(std::min(nxc,nyc) >= min_coarse);

		return ng;
	}
	void setup_space(int nlevels)
	{
		auto params = build_kernel_params(*conf);

		for (auto i : range(nlevels-1)) {
			auto & fop = levels.get<i>.A;
			auto nx = fop.shape(0);
			auto ny = fop.shape(1);
			auto nxc = (nx-1)/2. + 1;
			auto nyc = (ny-1)/2. + 1;

			levels.add(nxc, nyc);

			log::debug << "Created coarse grid with dimensions: " << nxc
			           << ", " << nyc << std::endl;
		}

		auto & cop = levels.get<levels.size()-1>().A;
		auto nxc = cop.shape(0);
		auto nyc = cop.shape(1);

		len_t abd_len_0 = nxc+2;
		if (params->periodic[0] or params->periodic[1])
			abd_len_0 = nxc*nyc;
		ABD = grid_func(abd_len_0, nxc*nyc, 0);
		bbd = new real_t[ABD.len(1)];
	}
};

}}

#endif
