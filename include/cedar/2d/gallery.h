#ifndef CEDAR_2D_GALLERY_H
#define CEDAR_2D_GALLERY_H

#include <cedar/types.h>
#include <cedar/2d/stencil_op.h>

namespace cedar { namespace cdr2 { namespace gallery {

stencil_op poisson(cedar::len_t nx, cedar::len_t ny);

stencil_op diag_diffusion(cedar::len_t nx, cedar::len_t ny,
                                 cedar::real_t dx, cedar::real_t dy);

stencil_op fe(cedar::len_t nx, cedar::len_t ny);

}}}

#endif
