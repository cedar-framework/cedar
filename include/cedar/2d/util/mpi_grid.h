#ifndef CEDAR_2D_UTIL_MPI_H
#define CEDAR_2D_UTIL_MPI_H

#include <mpi.h>
#include "cedar/types.h"
#include "cedar/2d/types.h"
#include "cedar/mpi/grid_topo.h"


namespace cedar { namespace cdr2 { namespace util { namespace mpi {

bool has_boundary(grid_topo & grid, cdr2::dir dir);

}}}}
#endif
