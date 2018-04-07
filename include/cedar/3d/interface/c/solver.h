#ifndef CEDAR_3D_INTERFACE_SOLVER_H
#define CEDAR_3D_INTERFACE_SOLVER_H

#include <cedar/3d/interface/c/operator.h>

#ifdef __cplusplus
extern "C" {
#endif

struct bmg3_slv;
typedef struct bmg3_slv* bmg3_solver;

bmg3_solver bmg3_solver_create(bmg3_operator *op);
void bmg3_solver_run(bmg3_solver op, double *x, const double *b);
void bmg3_solver_destroy(bmg3_solver);


#ifdef __cplusplus
}
#endif

#endif
