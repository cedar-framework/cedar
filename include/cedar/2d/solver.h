#ifndef CEDAR_2D_CORE_SOLVER_H
#define CEDAR_2D_CORE_SOLVER_H

#include <array>
#include <functional>
#include <cstddef>
#include <algorithm>

#include <cedar/level.h>
#include <cedar/multilevel.h>
#include <cedar/config.h>
#include <cedar/2d/level_container.h>
#include <cedar/2d/stencil_op.h>
#include <cedar/2d/relax_stencil.h>
#include <cedar/2d/prolong_op.h>
#include <cedar/2d/restrict_op.h>
#include <cedar/2d/kernel_manager.h>
#include <cedar/services/mempool.h>

namespace cedar { namespace cdr2 {

template<class sten>
struct level2 : public level<sten, stypes>
{
	using parent = level<sten, stypes>;
	level2(services::mempool & mpool, len_t nx, len_t ny) : parent::level(nx, ny)
	{
		this->x = grid_func(nx, ny);
		this->res = grid_func(nx, ny);
		this->res.set(0.0);
		this->b = grid_func(nx, ny);
		this->SOR[0] = relax_stencil(nx, ny);
		this->SOR[1] = relax_stencil(nx, ny);
		this->R.associate(&this->P);
	}
	level2(services::mempool & mpool, stencil_op<sten> & A) : parent::level(A)
	{
		this->res = grid_func(A.shape(0), A.shape(1));
		this->res.set(0.0);
		this->SOR[0] = relax_stencil(A.shape(0), A.shape(1));
		this->SOR[1] = relax_stencil(A.shape(0), A.shape(1));
	}
};

template<class fsten>
class solver: public multilevel<exec_mode::serial, level_container<level2,fsten>, fsten, solver<fsten>>
{
public:
	using parent = multilevel<exec_mode::serial, level_container<level2, fsten>, fsten, solver<fsten>>;
	using parent::levels;
solver(stencil_op<fsten> & fop) : parent::multilevel(fop)
	{
		this->kman = build_kernel_manager(*this->conf);
		setup(fop);
	}
	solver(stencil_op<fsten> & fop,
	       std::shared_ptr<config> conf) :
	parent::multilevel(fop, conf)
	{
		this->kman = build_kernel_manager(*this->conf);
		setup(fop);
	}
	~solver() { delete[] this->bbd; }

	void setup(stencil_op<fsten> & fop)
	{
		parent::setup(fop);
		auto kparams = this->kman->get_params();
		if (kparams->offload)
			prefetch();
	}

	void prefetch()
	{
		auto & lvl = levels.template get<fsten>(0);
		memory::prefetch(lvl.A.data(), lvl.A.size());
		memory::prefetch(lvl.res.data(), lvl.res.size());
		memory::prefetch(lvl.SOR[0].data(), lvl.SOR[0].size());
		memory::prefetch(lvl.SOR[1].data(), lvl.SOR[1].size());

		for (std::size_t i = 1; i < levels.size() - 1; i++) {
			auto & lvl = levels.get(i);
			memory::prefetch(lvl.A.data(), lvl.A.size());
			memory::prefetch(lvl.res.data(), lvl.res.size());
			memory::prefetch(lvl.SOR[0].data(), lvl.SOR[0].size());
			memory::prefetch(lvl.SOR[1].data(), lvl.SOR[1].size());
			memory::prefetch(lvl.x.data(), lvl.x.size());
			memory::prefetch(lvl.b.data(), lvl.b.size());
			memory::prefetch(lvl.P.data(), lvl.P.size());
		}

		// memory::prefetch(this->ABD.data(), this->ABD.size());
		// memory::prefetch(this->bbd, this->ABD.len(1));
		memory::sync();
		log::status << "Prefetched: " << memory::prefetchsize() / 1e9 << " GB" << std::endl;
	}


	std::size_t compute_num_levels(stencil_op<fsten> & fop)
	{
		float nxc, nyc;
		int ng = 0;
		auto min_coarse = this->settings.min_coarse;

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
		service_manager<stypes> & sman = this->kman->services();
		this->levels.init(sman.get<mempool>(), nlevels);
		auto params = build_kernel_params(*(this->conf));

		for (auto i : range<std::size_t>(nlevels-1)) {
			// TODO: remove copy-paste coding
			if (i == 0) {
				auto & fop = this->levels.template get<fsten>(i).A;
				auto nx = fop.shape(0);
				auto ny = fop.shape(1);
				len_t nxc = (nx-1)/2. + 1;
				len_t nyc = (ny-1)/2. + 1;

				this->levels.add(sman.get<mempool>(), nxc, nyc);

				log::debug << "Created coarse grid with dimensions: " << nxc
				           << ", " << nyc << std::endl;
			} else {
				auto & fop = this->levels.get(i).A;
				auto nx = fop.shape(0);
				auto ny = fop.shape(1);
				len_t nxc = (nx-1)/2. + 1;
				len_t nyc = (ny-1)/2. + 1;

				this->levels.add(sman.get<mempool>(), nxc, nyc);

				log::debug << "Created coarse grid with dimensions: " << nxc
				           << ", " << nyc << std::endl;
			}
		}

		auto & cop = this->levels.get(this->levels.size()-1).A;
		auto nxc = cop.shape(0);
		auto nyc = cop.shape(1);

		len_t abd_len_0 = nxc+2;
		if (params->periodic[0] or params->periodic[1] or params->offload)
			abd_len_0 = nxc*nyc;
		this->ABD = grid_func(abd_len_0, nxc*nyc, 0);
		this->bbd = new real_t[this->ABD.len(1)];
	}

	void give_op(std::unique_ptr<stencil_op<fsten>> fop) {fop_ref = std::move(fop);}

protected:
	std::unique_ptr<stencil_op<fsten>> fop_ref;
};

}}

#endif
