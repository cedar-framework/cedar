import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt('sol.txt')

plt.imshow(a)
plt.colorbar()
plt.show()
