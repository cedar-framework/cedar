#ifndef BOXMG_2D_UTIL_TOPO_H
#define BOXMG_2D_UTIL_TOPO_H

#include <mpi.h>
#include "boxmg-common.h"
#include "core/mpi/grid_topo.h"

namespace boxmg { namespace bmg2d { namespace util {

			bmg2d::mpi::topo_ptr create_topo(MPI_Comm comm, len_t nx, len_t ny);

}}}
#endif
