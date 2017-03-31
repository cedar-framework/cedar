#ifndef CEDAR_2D_INTERFACE_TOPO_H
#define CEDAR_2D_INTERFACE_TOPO_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct bmg2_topology;
typedef struct bmg2_topology* bmg2_topo;

bmg2_topo bmg2_topo_create(MPI_Comm comm,
                           unsigned int ngx,
                           unsigned int ngy,
                           unsigned int lnx[],
                           unsigned int lny[],
                           int nprocx,
                           int nprocy);


#ifdef __cplusplus
}
#endif

#endif
