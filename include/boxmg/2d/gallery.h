#ifndef BOXMG_2D_GALLERY_H
#define BOXMG_2D_GALLERY_H

#include <boxmg/types.h>
#include <boxmg/2d/stencil_op.h>

namespace boxmg { namespace bmg2d {

stencil_op create_poisson(boxmg::len_t nx, boxmg::len_t ny);

}}

#endif
