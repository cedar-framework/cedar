#ifndef CEDAR_3D_MPI_GALLERY_H
#define CEDAR_3D_MPI_GALLERY_H

#include <cedar/types.h>
#include <cedar/3d/mpi/stencil_op.h>

namespace cedar { namespace cdr3 { namespace mpi { namespace gallery {

stencil_op<seven_pt> poisson(topo_ptr grid);
stencil_op<xxvii_pt> fe(topo_ptr grid);

}}}}


#endif
