#ifndef CEDAR_2D_MPI_GALLERY_H
#define CEDAR_2D_MPI_GALLERY_H

#include <cedar/2d/mpi/stencil_op.h>

namespace cedar { namespace cdr2 { namespace mpi { namespace gallery {

stencil_op<five_pt> poisson(topo_ptr grid);
stencil_op<five_pt> diag_diffusion(topo_ptr grid, real_t dx, real_t dy);
stencil_op<nine_pt> fe(topo_ptr grid);

}}}}

#endif
