#ifndef CEDAR_UTIL_IO_H
#define CEDAR_UTIL_IO_H

#include <fstream>

#include <cedar/3d/mpi/stencil_op.h>

namespace cedar { namespace cdr3 { namespace util {

template<class sten>
void writeascii(mpi::stencil_op<sten>  & so);

template<class sten>
void readascii(mpi::stencil_op<sten> & so);

}}}
#endif
