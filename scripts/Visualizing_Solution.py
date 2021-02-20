import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

T = np.genfromtxt('scripts/T_solution.csv', delimiter=',')

N, M = T.shape

fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
t = np.linspace(0, N*1E-3, N+1)
x = np.linspace(0, 6.35E-3, M+1)

X, Y = np.meshgrid(x, y, indexing='ij')

# Plot the surface.
surf = ax.plot_surface(X, Y, T, cmap='magma')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()