#ifndef BOXMG_2D_GALLERY_H
#define BOXMG_2D_GALLERY_H

#include <boxmg/types.h>
#include <boxmg/2d/stencil_op.h>

namespace boxmg { namespace bmg2d { namespace gallery {

stencil_op poisson(boxmg::len_t nx, boxmg::len_t ny);

stencil_op diag_diffusion(boxmg::len_t nx, boxmg::len_t ny,
                                 boxmg::real_t dx, boxmg::real_t dy);

stencil_op fe(boxmg::len_t nx, boxmg::len_t ny);

}}}

#endif
