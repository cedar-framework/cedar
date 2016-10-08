#ifndef BOXMG_3D_INTERFACE_OPERATOR_H
#define BOXMG_3D_INTERFACE_OPERATOR_H

#include <boxmg/3d/interface/c/topo.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <boxmg/3d/base_types.h>

struct bmg3_op;
typedef struct bmg3_op* bmg3_operator;

typedef struct {
	unsigned int i;
	unsigned int j;
	unsigned int k;
	bmg3_dir dir;
} grid_coord_3d;

bmg3_operator bmg3_operator_create(bmg3_topo topo);

bmg3_operator bmg3_operator_create_sevenpt(bmg3_topo topo);

void bmg3_operator_set(bmg3_operator, unsigned int nvals, grid_coord_3d coords[], double vals[]);

void bmg3_operator_apply(bmg3_operator, const double *x, double *b);

void bmg3_operator_dump(bmg3_operator);

void bmg3_operator_destroy(bmg3_operator);

#ifdef __cplusplus
}
#endif

#endif
