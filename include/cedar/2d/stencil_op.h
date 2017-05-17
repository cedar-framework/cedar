#ifndef CEDAR_2D_STENCIL_OP_H
#define CEDAR_2D_STENCIL_OP_H

#include <cedar/stencil_op_nd.h>
#include <cedar/2d/base_types.h>

namespace cedar { namespace cdr2 {

		enum class five_pt {c=BMG2_C, w=BMG2_W, s=BMG2_S, ndirs};
		enum class nine_pt {c=BMG2_C, w=BMG2_W, s=BMG2_S, sw=BMG2_SW, nw=BMG2_NW, ndirs};

		template<class stype>
			using stencil_op = stencil_op_nd<2, stype>;

}}

#endif
