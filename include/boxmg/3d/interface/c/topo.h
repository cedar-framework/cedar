#ifndef BOXMG_3D_INTERFACE_TOPO_H
#define BOXMG_3D_INTERFACE_TOPO_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct bmg3_topology;
typedef struct bmg3_topology* bmg3_topo;

bmg3_topo bmg3_topo_create(MPI_Comm comm,
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
