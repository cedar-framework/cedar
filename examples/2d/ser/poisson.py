from cedar.cdr2 import *
import numpy as np
import matplotlib.pyplot as plt

def create_poisson(nx, ny):
    so = stencil_op_five(nx, ny)
    a = so.toarray()

    hx = 1.0 / (nx + 1)
    hy = 1.0 / (ny + 1)
    xh = hy / hx
    yh = hx / hy
    i1 = nx + 1
    j1 = ny + 1

    a[1:i1, 2:j1, five_pt.s] = yh
    a[2:i1, 1:j1, five_pt.w] = xh
    a[1:i1, 1:j1, five_pt.c] = 2*xh + 2*yh

    return so


def rhs(x, y):
    return 8*(np.pi*np.pi)*np.sin(2*np.pi*x)*np.sin(2*np.pi*y)


nx = 101
ny = nx

x = np.linspace(0, 1, (nx+2))
y = np.linspace(0, 1, (ny+2))
xv, yv = np.meshgrid(x, y)


hx = 1.0 / (nx + 1)
hy = 1.0 / (ny + 1)
h2 = hx * hy

b = grid_func(nx, ny)
x = grid_func(nx, ny)

b.toarray()[:,:] = h2 * rhs(xv, yv)

#a = gallery.poisson(nx, ny)
a = create_poisson(nx, ny)

solve(a, x, b)

plt.pcolormesh(xv, yv, x.toarray())

plt.show()
