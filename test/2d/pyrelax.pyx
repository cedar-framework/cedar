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

    inds = np.zeros(nx*ny, dtype=np.int32)
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


cdef public double get_high_freq(int nx, int ny, double *src):
    x = np.zeros((ny, nx))
    for j in range(1, ny+1):
        for i in range(1, nx+1):
            x[j-1][i-1] = src[j*(nx+2) + i]

    # import matplotlib as mpl
    # mpl.use('Qt5Agg')
    # from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    # X = np.linspace(0, 1, nx)
    # Y = np.linspace(0, 1, ny)
    # X,Y = np.meshgrid(X, Y)

    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # surf = ax.plot_wireframe(X, Y, x)
    # plt.show()
    freq0 = np.fft.fft(x[ny//4, :])
    freq1 = np.fft.fft(x[:, nx//4])
    switch = 2
    mx = np.max(np.abs(freq0[nx//(switch*2):nx//2]))
    my = np.max(np.abs(freq1[ny//(switch*2):ny//2]))
    return max(mx,my)
    # plt.plot(np.abs(freq0[:nx//2]), label='0')
    # plt.plot(np.abs(freq1[:ny//2]), label='1')
    # plt.legend(loc=0)
    # plt.show()
