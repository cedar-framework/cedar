from libc.stdlib cimport malloc
import numpy as np
import pyamg
import scipy.sparse as sparse
import scipy.sparse.linalg as sla
cimport numpy as np

cdef public double gs_iter(int nx, int ny, int nsweeps):
    A = pyamg.gallery.poisson((nx,ny))
    b = np.zeros(nx*ny)
    x = np.random.rand(nx*ny)

    r = b - A.dot(x)
    n = x.shape[0]
    D = A.diagonal()[0]*sparse.eye(n, format='csr')
    L = sparse.tril(A, -1, format='csr')
    U = sparse.triu(A, 1, format='csr')
    DLInv = sla.inv(D+L)
    G = DLInv * U
    for i in range(nsweeps):
        x = (DLInv*b) - (G*x)

    r = b - A.dot(x)
    return np.linalg.norm(r, np.inf)
