from pycedar import *
import numpy as np

def stencil_op_to_array(self):
    return np.array(self, copy=False)

for opclass in [stencil_op_five, stencil_op_nine]:
    opclass.toarray = stencil_op_to_array

