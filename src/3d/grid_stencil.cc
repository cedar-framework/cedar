#include <boxmg/types.h>
#include <boxmg/util/log.h>
#include <boxmg/3d/grid_stencil.h>

using namespace boxmg::bmg3;

grid_stencil::grid_stencil(len_t nx, len_t ny, len_t nz, unsigned int nghosts, bool intergrid, bool symmetric, bool five_pt) : array(), symmetric(not intergrid and symmetric), five_pt_(five_pt)
{
	num_ghosts = nghosts;
	std::vector<len_t> lens = {nx, ny, nz};
	for (auto i : boxmg::range(3)) {
		range_[i] = boxmg::range(static_cast<len_t>(nghosts), static_cast<len_t>(lens[i]+nghosts));
		grange_[i] = boxmg::range(static_cast<len_t>(0), lens[i] + 2*nghosts);
		lens[i] += nghosts*2;
	}

	if (this->symmetric) {
		reshape(lens[0], lens[1], lens[2], 14);
	} else {
		if (intergrid)
			reshape(lens[0], lens[1], lens[2], 26);
		else
			log::error << "non-symmetric problems not supported" << std::endl;
	}
}


grid_stencil & grid_stencil::operator=(grid_stencil&& gs)
{
	vec = std::move(gs.vec);
	strides = std::move(gs.strides);
	extents = std::move(gs.extents);
	five_pt_ = gs.five_pt_;
	symmetric = gs.symmetric;
	num_ghosts = gs.num_ghosts;
	range_ = std::move(gs.range_);
	grange_ = std::move(gs.grange_);

	return *this;
}
