from _cdr2 import *
import gallery
import numpy as np

def bufobj_to_array(self):
    return np.array(self, copy=False)

for bufobj in [stencil_op_five, stencil_op_nine, grid_func]:
    bufobj.toarray = bufobj_to_array

del bufobj
del bufobj_to_array
