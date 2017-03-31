#ifndef CEDAR_2D_INTERFACE_SOLVER_H
#define CEDAR_2D_INTERFACE_SOLVER_H

#include <cedar/2d/interface/c/operator.h>

#ifdef __cplusplus
extern "C" {
#endif

struct bmg2_slv;
typedef struct bmg2_slv* bmg2_solver;

bmg2_solver bmg2_solver_create(bmg2_operator *op);
void bmg2_solver_run(bmg2_solver op, double *x, const double *b);
void bmg2_solver_destroy(bmg2_solver);


#ifdef __cplusplus
}
#endif

#endif
