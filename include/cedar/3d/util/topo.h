#ifndef CEDAR_3D_UTIL_TOPO_H
#define CEDAR_3D_UTIL_TOPO_H

#include <mpi.h>
#include <cedar/types.h>
#include <cedar/mpi/grid_topo.h>
#include <cedar/config.h>

namespace cedar { namespace cdr3 { namespace util {

			topo_ptr create_topo(MPI_Comm comm, len_t nx, len_t ny, len_t nz);
			topo_ptr create_topo(MPI_Comm comm, len_t npx, len_t npy, len_t npz,
			                     len_t nlx, len_t nly, len_t nlz);
			topo_ptr create_topo_global(MPI_Comm comm, len_t ngx, len_t ngy, len_t ngz);
			topo_ptr create_topo(config & conf);
			topo_ptr model_topo(int np, len_t ngx, len_t ngy, len_t ngz);
			topo_ptr coarsen_topo(topo_ptr topof);
}}}

#endif
