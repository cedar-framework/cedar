#ifndef BOXMG_2D_MPI_GALLERY_H
#define BOXMG_2D_MPI_GALLERY_H

#include <boxmg/2d/mpi/stencil_op.h>

namespace boxmg { namespace bmg2d { namespace mpi { namespace gallery {

stencil_op poisson(topo_ptr grid);

}}}}

#endif
