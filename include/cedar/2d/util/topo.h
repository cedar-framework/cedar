#ifndef CEDAR_2D_UTIL_TOPO_H
#define CEDAR_2D_UTIL_TOPO_H

#include <mpi.h>
#include <cedar/types.h>
#include <cedar/config/reader.h>
#include <cedar/mpi/grid_topo.h>

namespace cedar { namespace cdr2 { namespace util {
			topo_ptr create_topo(MPI_Comm comm, len_t nx, len_t ny);
			topo_ptr create_topo(MPI_Comm comm, int npx, int npy, len_t nlx, len_t nly);
			topo_ptr create_topo(config::reader & conf);
			topo_ptr model_topo(int nprocx, int nprocy, len_t ngx, len_t ngy);
			topo_ptr coarsen_topo(topo_ptr topof);
			topo_ptr create_topo_global(MPI_Comm comm, len_t ngx, len_t ngy);
}}}
#endif
