from libc.stdlib cimport malloc
import numpy as np
import pyamg
from pyamg.relaxation.relaxation import gauss_seidel_indexed
import scipy.sparse as sparse
import scipy.sparse.linalg as sla
cimport numpy as np

cdef public void gs_iter(int nx, int ny, int nsweeps, int npt_stencil, double *dst):
    type = 'FD'
    if npt_stencil == 9:
        type = 'FE'

    A = pyamg.gallery.poisson((nx,ny), type=type)
    b = np.zeros(nx*ny)
    x = np.ones(nx*ny)

    inds = np.zeros(nx*ny, dtype=np.int)
    cnt = 0
    if npt_stencil == 9:
        for j in range(0, ny, 2):
            for i in range(0, nx, 2):
                inds[cnt] = j*nx + i
                cnt += 1
            for i in range(1, nx, 2):
                inds[cnt] = j*nx + i
                cnt += 1

        for j in range(1, ny, 2):
            for i in range(0, nx, 2):
                inds[cnt] = j*nx + i
                cnt += 1
            for i in range(1, nx, 2):
                inds[cnt] = j*nx + i
                cnt += 1
    else:
        for k in range(2):
            for j in range(0, ny):
                if j % 2 == k:
                    for i in range(0, nx, 2):
                        inds[cnt] = j*nx+i
                        cnt+=1
                else:
                    for i in range(1, nx, 2):
                        inds[cnt] = j*nx+i
                        cnt+=1


    gauss_seidel_indexed(A, x, b, inds, iterations=nsweeps, sweep='forward')
    gauss_seidel_indexed(A, x, b, inds, iterations=nsweeps, sweep='backward')

    for i in range(nx*ny):
        dst[i] = x[i]
