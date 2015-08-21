#ifndef BOXMG_2D_INTERFACE_SOLVER_H
#define BOXMG_2D_INTERFACE_SOLVER_H


struct bmg2_instance;
typedef struct bmg2_instance bmg2_instance;

struct bmg2_solver
{
	bmg2_instance * slv;
};

typedef struct bmg2_solver * bmg2_solver;

struct bmg2_solver;
typedef struct bmg2_solver bmg2_solver;

struct bmg2_op;
typedef struct bmg2_op* bmg2_op;

bmg2_solver bmg2_solver_create(void);
void bmg2_solver_destroy(bmg2_solver);

#endif
