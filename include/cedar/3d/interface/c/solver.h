#ifndef CEDAR_3D_INTERFACE_SOLVER_H
#define CEDAR_3D_INTERFACE_SOLVER_H

#include <cedar/3d/interface/c/operator.h>

#ifdef __cplusplus
extern "C" {
#endif

struct cdr3_slv;
typedef struct cdr3_slv* cdr3_solver;

cdr3_solver cdr3_solver_create(cdr3_operator *op);
void cdr3_solver_run(cdr3_solver op, double *x, const double *b);
void cdr3_solver_destroy(cdr3_solver);


#ifdef __cplusplus
}
#endif

#endif
