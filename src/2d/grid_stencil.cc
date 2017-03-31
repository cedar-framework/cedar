#include <stdlib.h>

#include <cedar/types.h>
#include <cedar/2d/grid_stencil.h>

using namespace cedar::cdr2;


grid_stencil::grid_stencil(len_t nx, len_t ny, unsigned int nghosts, bool intergrid, bool symmetric, bool five_pt): array(), symmetric(not intergrid and symmetric), five_pt_(five_pt)
{
	num_ghosts = nghosts;
	range_[0] = cedar::range(static_cast<len_t>(nghosts), static_cast<len_t>(nx + nghosts));
	range_[1] = cedar::range(static_cast<len_t>(nghosts), static_cast<len_t>(ny + nghosts));

	grange_[0] = cedar::range(static_cast<len_t>(0), nx + 2*nghosts);
	grange_[1] = cedar::range(static_cast<len_t>(0), ny + 2*nghosts);

	nx += nghosts*2;
	ny += nghosts*2;

	if (this->symmetric) {
		reshape(nx, ny, 5);
		//size = (nx+2)*(ny+2)*5;  // temporarily add space for mpi stencil
	} else {
		reshape(nx,ny,8);
		//size = (nx + 1) * (ny + 1) * 8;
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
