#ifndef CEDAR_RESCUDA_H
#define CEDAR_RESCUDA_H


void residual_cuda9(int ilen, int jlen, double *r, const double *so, const double *x, const double *b);
void set_constant_cuda(double *arr, int arrsize, double val);

#endif
