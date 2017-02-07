#ifndef BOXMG_3D_GALLERY_H
#define BOXMG_3D_GALLERY_H

#include <boxmg/types.h>
#include <boxmg/3d/stencil_op.h>

namespace boxmg { namespace bmg3 { namespace gallery {

stencil_op poisson(boxmg::len_t nx, boxmg::len_t ny, boxmg::len_t nz);

}}}


#endif
