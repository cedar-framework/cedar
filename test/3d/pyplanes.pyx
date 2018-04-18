from libc.stdlib cimport malloc
import numpy as np
import pyamg
import scipy.sparse as sparse
import scipy.sparse.linalg as sla
cimport numpy as np

cdef public void pyplane(int nx, int ny, int nz, double *cx, double *cb, int ipl):
    b = np.zeros(nx * ny)

    for j in range(ny):
        for i in range(nx):
            b[j*nx+i] = cb[(j+1)*(nx+2) + (i+1)]

    hx = 1. / (nx + 1)
    hy = 1. / (ny + 1)
    hz = 1. / (nz + 1)
    xh = hy*hz/hx
    yh = hx*hz/hy
    zh = hx*hy/hz

    stencil = [[0,     -1*yh, 0],
               [-1*xh, 2*xh + 2*yh + 2*zh, -1*xh],
               [0, -1 * yh, 0]]
    A = pyamg.gallery.stencil_grid(stencil, (ny, nx), dtype='float', format='csr')
    x = sla.spsolve(A, b)

    for j in range(ny):
        for i in range(nx):
            idx = ipl * (nx+2)*(ny+2) + (j+1) * (nx+2) + (i+1)
            cx[idx] = x[j*nx+i]


cdef public void pyplane27(int nx, int ny, int nz, double *cx, double *cb, int ipl):
    b = np.zeros(nx * ny)

    for j in range(ny):
        for i in range(nx):
            b[j*nx+i] = cb[(j+1)*(nx+2) + (i+1)]

    stencil = [[-1.,-1.,-1.,],
               [-1., 26., -1.],
               [-1., -1., -1.]]
    A = pyamg.gallery.stencil_grid(stencil, (ny, nx), dtype='float', format='csr')
    x = sla.spsolve(A, b)

    for j in range(ny):
        for i in range(nx):
            idx = ipl * (nx+2)*(ny+2) + (j+1) * (nx+2) + (i+1)
            cx[idx] = x[j*nx+i]


cdef public void pyplane_xz(int nx, int ny, int nz, double *cx, double *cb, int ipl):
    b = np.zeros(nx * nz)

    for k in range(nz):
        for i in range(nx):
            b[k*nx+i] = cb[(k+1)*(nx+2) + (i+1)]

    hx = 1. / (nx + 1)
    hy = 1. / (ny + 1)
    hz = 1. / (nz + 1)
    xh = hy*hz/hx
    yh = hx*hz/hy
    zh = hx*hy/hz

    stencil = [[0,     -1*zh, 0],
               [-1*xh, 2*xh + 2*yh + 2*zh, -1*xh],
               [0, -1 * zh, 0]]
    A = pyamg.gallery.stencil_grid(stencil, (nz, nx), dtype='float', format='csr')
    x = sla.spsolve(A, b)

    for k in range(nz):
        for i in range(nx):
            idx = (k+1) * (nx+2)*(ny+2) + ipl * (nx+2) + (i+1)
            cx[idx] = x[k*nx+i]


cdef public void pyplane_xz27(int nx, int ny, int nz, double *cx, double *cb, int ipl):
    b = np.zeros(nx * nz)

    for k in range(nz):
        for i in range(nx):
            b[k*nx+i] = cb[(k+1)*(nx+2) + (i+1)]

    stencil = [[-1.,-1.,-1.,],
               [-1., 26., -1.],
               [-1., -1., -1.]]

    A = pyamg.gallery.stencil_grid(stencil, (nz, nx), dtype='float', format='csr')
    x = sla.spsolve(A, b)

    for k in range(nz):
        for i in range(nx):
            idx = (k+1) * (nx+2)*(ny+2) + ipl * (nx+2) + (i+1)
            cx[idx] = x[k*nx+i]


cdef public void pyplane_yz(int nx, int ny, int nz, double *cx, double *cb, int ipl):
    b = np.zeros(ny * nz)

    for k in range(nz):
        for j in range(ny):
            b[k*ny+j] = cb[(k+1)*(ny+2) + (j+1)]

    hx = 1. / (nx + 1)
    hy = 1. / (ny + 1)
    hz = 1. / (nz + 1)
    xh = hy*hz/hx
    yh = hx*hz/hy
    zh = hx*hy/hz

    stencil = [[0,     -1*zh, 0],
               [-1*yh, 2*xh + 2*yh + 2*zh, -1*yh],
               [0, -1 * zh, 0]]
    A = pyamg.gallery.stencil_grid(stencil, (nz, ny), dtype='float', format='csr')
    x = sla.spsolve(A, b)

    for k in range(nz):
        for j in range(ny):
            idx = (k+1) * (nx+2)*(ny+2) + (j+1) * (nx+2) + ipl
            cx[idx] = x[k*ny+j]


cdef public void pyplane_yz27(int nx, int ny, int nz, double *cx, double *cb, int ipl):
    b = np.zeros(ny * nz)

    for k in range(nz):
        for j in range(ny):
            b[k*ny+j] = cb[(k+1)*(ny+2) + (j+1)]

    stencil = [[-1.,-1.,-1.,],
               [-1., 26., -1.],
               [-1., -1., -1.]]
    A = pyamg.gallery.stencil_grid(stencil, (nz, ny), dtype='float', format='csr')
    x = sla.spsolve(A, b)

    for k in range(nz):
        for j in range(ny):
            idx = (k+1) * (nx+2)*(ny+2) + (j+1) * (nx+2) + ipl
            cx[idx] = x[k*ny+j]
