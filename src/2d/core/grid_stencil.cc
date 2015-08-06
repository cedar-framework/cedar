#include <stdlib.h>

#include "grid_stencil.h"

using namespace boxmg::bmg2d::core;


GridStencil::GridStencil(len_t nx, len_t ny, unsigned int nghosts, bool intergrid, bool symmetric, bool five_pt): Array(nx,ny,nghosts,false), symmetric(not intergrid and symmetric), five_pt_(five_pt)
{
	len_t size;
	nx += nghosts*2;
	ny += nghosts*2;
	if (this->symmetric) {
		size = (nx+1)*(ny+1)*5;
		//size = (nx+2)*(ny+2)*5;  // temporarily add space for mpi stencil
	} else {
		size = nx*ny*8;
		//size = (nx + 1) * (ny + 1) * 8;
	}

	data_.resize(size);

	ostride = len_[0]*len_[1];
}
