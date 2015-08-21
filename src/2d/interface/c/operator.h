#ifndef BOXMG_2D_INTERFACE_OPERATOR_H
#define BOXMG_2D_INTERFACE_OPERATOR_H

#ifdef __cplusplus
extern "C" {
#endif

struct bmg2_op;
typedef struct bmg2_op* bmg2_operator;

bmg2_operator bmg2_operator_create(unsigned int nx, unsigned int ny);
void bmg2_operator_test(bmg2_operator);
void bmg2_operator_destroy(bmg2_operator);

#ifdef __cplusplus
}
#endif

#endif
