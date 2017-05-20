#ifndef CEDAR_2D_CORE_SOLVER_H
#define CEDAR_2D_CORE_SOLVER_H

#include <array>
#include <functional>
#include <cstddef>
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
level_container(stencil_op<sten> & fine_op) : fine_op(fine_op) {}
	void init(std::size_t nlevels);
	void add(len_t nx, len_t ny) {
		lvls_nine.emplace_back(nx, ny);
	}
	template<class rsten=nine_pt> level<rsten>&  get(std::size_t i);
	std::size_t size() { return lvls_nine.size() + lvls_five.size(); }

protected:
	stencil_op<sten> & fine_op;
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

template<> void level_container<nine_pt>::init(std::size_t nlevels)
{
	lvls_nine.reserve(nlevels);
	lvls_nine.emplace_back(fine_op);
}

template<> void level_container<five_pt>::init(std::size_t nlevels)
{
	lvls_five.emplace_back(fine_op);
	lvls_nine.reserve(nlevels-1);
}

template<class fsten>
class solver: public multilevel<level, level_container<fsten>,
	stencil_op, grid_func, kernel::registry, fsten, solver<fsten>>
{
public:
	using parent = multilevel<level, level_container<fsten>,
		stencil_op, grid_func, kernel::registry, fsten, solver<fsten>>;
solver(stencil_op<fsten> & fop) : parent::multilevel(fop)
	{
		this->kreg = std::make_shared<kernel::registry>(*(this->conf));
		parent::setup(fop);
	}
	solver(stencil_op<fsten> & fop,
	       std::shared_ptr<config::reader> conf) :
	parent::multilevel(fop, conf)
	{
		this->kreg = std::make_shared<kernel::registry>(*(this->conf));
		parent::setup(fop);
	}
	~solver() { delete[] this->bbd; }
	std::size_t compute_num_levels(stencil_op<fsten> & fop)
	{
		float nxc, nyc;
		int ng = 0;
		auto min_coarse = this->conf->template get<std::size_t>("solver.min-coarse", 3);

		auto nx = fop.shape(0);
		auto ny = fop.shape(1);

		do {
			ng++;
			nxc = (nx-1)/(1<<ng) + 1;
			nyc = (ny-1)/(1<<ng) + 1;
		} while(std::min(nxc,nyc) >= min_coarse);

		return ng;
	}

	void setup_space(std::size_t nlevels)
	{
		this->levels.init(nlevels);
		auto params = build_kernel_params(*(this->conf));

		for (auto i : range<std::size_t>(nlevels-1)) {
			// TODO: remove copy-paste coding
			if (i == 0) {
				auto & fop = this->levels.template get<fsten>(i).A;
				auto nx = fop.shape(0);
				auto ny = fop.shape(1);
				len_t nxc = (nx-1)/2. + 1;
				len_t nyc = (ny-1)/2. + 1;

				this->levels.add(nxc, nyc);

				log::debug << "Created coarse grid with dimensions: " << nxc
				           << ", " << nyc << std::endl;
			} else {
				auto & fop = this->levels.get(i).A;
				auto nx = fop.shape(0);
				auto ny = fop.shape(1);
				len_t nxc = (nx-1)/2. + 1;
				len_t nyc = (ny-1)/2. + 1;

				this->levels.add(nxc, nyc);

				log::debug << "Created coarse grid with dimensions: " << nxc
				           << ", " << nyc << std::endl;
			}
		}

		auto & cop = this->levels.get(this->levels.size()-1).A;
		auto nxc = cop.shape(0);
		auto nyc = cop.shape(1);

		len_t abd_len_0 = nxc+2;
		if (params->periodic[0] or params->periodic[1])
			abd_len_0 = nxc*nyc;
		this->ABD = grid_func(abd_len_0, nxc*nyc, 0);
		this->bbd = new real_t[this->ABD.len(1)];
	}
};

}}

#endif
