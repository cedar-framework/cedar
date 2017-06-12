#ifndef CEDAR_2D_CORE_SOLVER_H
#define CEDAR_2D_CORE_SOLVER_H

#include <array>
#include <functional>
#include <cstddef>
#include <algorithm>

#include <cedar/level.h>
#include <cedar/multilevel.h>
#include <cedar/config/reader.h>
#include <cedar/2d/level_container.h>
#include <cedar/2d/stencil_op.h>
#include <cedar/2d/relax_stencil.h>
#include <cedar/2d/inter/prolong_op.h>
#include <cedar/2d/inter/restrict_op.h>
#include <cedar/2d/kernel/registry.h>

namespace cedar { namespace cdr2 {

template<class sten>
struct level2 : public level<sten, stypes>
{
	using parent = level<sten, stypes>;
level2(len_t nx, len_t ny) : parent::level(nx, ny)
	{
		this->SOR = {{relax_stencil(nx, ny),
		              relax_stencil(nx, ny)}};
		this->R.associate(&this->P);
	}
level2(stencil_op<sten> & A) : parent::level(A)
	{
		this->res = grid_func(A.shape(0), A.shape(1));
		this->SOR = {{relax_stencil(A.shape(0), A.shape(1)),
		              relax_stencil(A.shape(0), A.shape(1))}};
	}
};

template<class fsten>
class solver: public multilevel<level_container<level2,fsten>,
	typename kernel::registry::parent, fsten, solver<fsten>>
{
public:
	using parent = multilevel<level_container<level2, fsten>,
		typename kernel::registry::parent, fsten, solver<fsten>>;
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

	void give_op(std::unique_ptr<stencil_op<fsten>> fop) {fop_ref = std::move(fop);}

protected:
	std::unique_ptr<stencil_op<fsten>> fop_ref;
};

template<class fsten>
	void solve(stencil_op<fsten> & A, grid_func & x, const grid_func & b)
{
	solver<fsten> bmg(A);
	bmg.solve(b, x);
}

}}

#endif
