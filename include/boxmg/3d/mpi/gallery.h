#ifndef BOXMG_3D_MPI_GALLERY_H
#define BOXMG_3D_MPI_GALLERY_H

#include <boxmg/types.h>
#include <boxmg/3d/mpi/stencil_op.h>

namespace boxmg { namespace bmg3 { namespace mpi { namespace gallery {

stencil_op poisson(topo_ptr grid);

}}}}


#endif
