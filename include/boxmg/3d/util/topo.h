#ifndef BOXMG_3D_UTIL_TOPO_H
#define BOXMG_3D_UTIL_TOPO_H

#include <mpi.h>
#include <boxmg/types.h>
#include <boxmg/mpi/grid_topo.h>

namespace boxmg { namespace bmg3 { namespace util {

			topo_ptr create_topo(MPI_Comm comm, len_t nx, len_t ny, len_t nz);

}}}

#endif
