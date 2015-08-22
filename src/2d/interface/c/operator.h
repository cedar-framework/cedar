#ifndef BOXMG_2D_INTERFACE_OPERATOR_H
#define BOXMG_2D_INTERFACE_OPERATOR_H

#ifdef __cplusplus
extern "C" {
#endif

#include "boxmg-2d/core/base_types.h"

struct bmg2_op;
typedef struct bmg2_op* bmg2_operator;

typedef struct {
	unsigned int i;
	unsigned int j;
	bmg2_dir dir;
} grid_coord;

bmg2_operator bmg2_operator_create(unsigned int nx, unsigned int ny);

void bmg2_operator_set(bmg2_operator, unsigned int nvals, grid_coord coords[], double vals[]);

void bmg2_operator_destroy(bmg2_operator);

#ifdef __cplusplus
}
#endif

#endif
