#ifndef CEDAR_3D_INTER_TYPES_H
#define CEDAR_3D_INTER_TYPES_H

namespace cedar { namespace cdr3 { namespace inter {

enum class dir {XYL=0, XYR=1, XYA=2, XYB=3, XZA=4, XZB=5,
		XYNE=6, XYSE=7, XYSW=8, XYNW=9, XZSW=10,
		XZNW=11,
		XZNE=12, XZSE=13, YZSW=14, YZNW=15, YZNE=16,
		YZSE=17, BSW=18, BNW=19, BNE=20, BSE=21,
		TSW=22, TNW=23, TNE=24, TSE=25, ndirs};

}}}
#endif
