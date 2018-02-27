#ifndef CEDAR_3D_STENCIL_OP_H
#define CEDAR_3D_STENCIL_OP_H

#include<iostream>

#include <cedar/types.h>
#include <cedar/stencil_op_nd.h>
#include <cedar/3d/base_types.h>

namespace cedar { namespace cdr3 {

		enum class seven_pt {p=BMG3_P, pw=BMG3_PW, ps=BMG3_PS, b=BMG3_B, ndirs};
		enum class xxvii_pt {p=BMG3_P, pw=BMG3_PW, ps=BMG3_PS, b=BMG3_B,
				psw=BMG3_PSW, pnw=BMG3_PNW, bw=BMG3_BW, bnw=BMG3_BNW,
				bn=BMG3_BN, bne=BMG3_BNE, be=BMG3_BE, bse=BMG3_BSE,
				bs=BMG3_BS, bsw=BMG3_BSW,ndirs};

		template<class stype>
			using stencil_op = stencil_op_nd<3, stype>;

}}
#endif
