import os
import sys

import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

T = np.genfromtxt(os.path.join(os.path.dirname(sys.path[0]), 'solutions/T_Solution.csv'), delimiter=',')

N, M = T.shape

fig = plt.figure()
ax = fig.gca(projection='3d')

t = np.linspace(0, (N-1)*1E-5, N)
x = np.linspace(0, 6.35E-3, M)

t_arr, x_arr = np.meshgrid(t, x, indexing='ij')

print(t_arr.shape, x_arr.shape, T.shape)

# Make data.
# Plot the surface.
surf = ax.plot_surface(t_arr, x_arr, T, cmap='magma')

ax.set_xlabel('t')
ax.set_ylabel('x')

plt.show()