#ifndef CEDAR_3D_INTERFACE_TOPO_H
#define CEDAR_3D_INTERFACE_TOPO_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct cdr3_topology;
typedef struct cdr3_topology* cdr3_topo;

cdr3_topo cdr3_topo_create(MPI_Comm comm,
                           unsigned int ngx,
                           unsigned int ngy,
                           unsigned int ngz,
                           unsigned int lnx[],
                           unsigned int lny[],
                           unsigned int lnz[],
                           int nprocx,
                           int nprocy,
                           int nprocz);


#ifdef __cplusplus
}
#endif

#endif
