#ifndef BOXMG_2D_UTIL_MPI_H
#define BOXMG_2D_UTIL_MPI_H

#include <mpi.h>
#include "boxmg-common.h"
#include "core/types.h"
#include "core/mpi/grid_topo.h"


namespace boxmg { namespace bmg2d { namespace util { namespace mpi {

bool has_boundary(bmg2d::mpi::grid_topo & grid, bmg2d::dir dir);

}}}}
#endif
