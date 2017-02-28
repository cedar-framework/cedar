#ifndef BOXMG_2D_INTERFACE_OPERATOR_H
#define BOXMG_2D_INTERFACE_OPERATOR_H

#include <boxmg/2d/interface/c/topo.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <boxmg/2d/base_types.h>

struct bmg2_op;
typedef struct bmg2_op* bmg2_operator;

typedef struct {
	unsigned int i;
	unsigned int j;
	bmg2_dir dir;
} grid_coord_2d;

bmg2_operator bmg2_operator_create(bmg2_topo topo);

bmg2_operator bmg2_operator_create_fivept(bmg2_topo topo);

void bmg2_operator_set(bmg2_operator, unsigned int nvals, grid_coord_2d coords[], double vals[]);

void bmg2_operator_apply(bmg2_operator, const double *x, double *b);

void bmg2_operator_dump(bmg2_operator);

void bmg2_operator_destroy(bmg2_operator);

#ifdef __cplusplus
}
#endif

#endif
