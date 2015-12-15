#include <stdlib.h>

#include <boxmg/types.h>
#include <boxmg/2d/grid_stencil.h>

using namespace boxmg::bmg2d;


grid_stencil::grid_stencil(len_t nx, len_t ny, unsigned int nghosts, bool intergrid, bool symmetric, bool five_pt): array(), symmetric(not intergrid and symmetric), five_pt_(five_pt), num_ghosts(nghosts)
{
	range_[0] = boxmg::range(static_cast<len_t>(nghosts), static_cast<len_t>(nx + nghosts));
	range_[1] = boxmg::range(static_cast<len_t>(nghosts), static_cast<len_t>(ny + nghosts));

	grange_[0] = boxmg::range(static_cast<len_t>(0), nx + 2*nghosts);
	grange_[1] = boxmg::range(static_cast<len_t>(0), ny + 2*nghosts);

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


const boxmg::range_t<boxmg::len_t> & grid_stencil::range(int i) const
{
		#ifdef DEBUG
		return range_.at(i);
		#else
		return range_[i];
		#endif
}


const boxmg::range_t<boxmg::len_t> & grid_stencil::grange(int i) const
{
	#ifdef DEBUG
	return grange_.at(i);
	#else
	return grange_[i];
	#endif
}


boxmg::len_t grid_stencil::shape(int i) const
{
	return this->len(i) - 2*num_ghosts;
}
