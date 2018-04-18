#ifndef CEDAR_3D_GALLERY_H
#define CEDAR_3D_GALLERY_H

#include <cedar/types.h>
#include <cedar/3d/stencil_op.h>

namespace cedar { namespace cdr3 { namespace gallery {

stencil_op<seven_pt> poisson(cedar::len_t nx, cedar::len_t ny, cedar::len_t nz);

stencil_op<seven_pt> diag_diffusion(cedar::len_t nx, cedar::len_t ny, cedar::len_t nz,
                                    cedar::real_t dx, cedar::real_t dy, cedar::real_t dz);
stencil_op<xxvii_pt> fe(cedar::len_t nx, cedar::len_t ny, cedar::len_t nz);

}}}


#endif
