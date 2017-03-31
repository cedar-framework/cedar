#ifndef CEDAR_3D_INTERFACE_OPERATOR_H
#define CEDAR_3D_INTERFACE_OPERATOR_H

#include <cedar/3d/interface/c/topo.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <cedar/3d/base_types.h>

struct cdr3_op;
typedef struct cdr3_op* cdr3_operator;

typedef struct {
	unsigned int i;
	unsigned int j;
	unsigned int k;
	cdr3_dir dir;
} grid_coord_3d;

cdr3_operator cdr3_operator_create(cdr3_topo topo);

cdr3_operator cdr3_operator_create_sevenpt(cdr3_topo topo);

void cdr3_operator_set(cdr3_operator, unsigned int nvals, grid_coord_3d coords[], double vals[]);

void cdr3_operator_apply(cdr3_operator, const double *x, double *b);

void cdr3_operator_dump(cdr3_operator);

void cdr3_operator_destroy(cdr3_operator);

#ifdef __cplusplus
}
#endif

#endif
