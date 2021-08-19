import numpy as np
import matplotlib.pyplot as plt

data = np.array(
    [
        [0.4,	0.9],
        [0.5,	1.2],
        [0.6,	1.3],
        [0.7,	1.4],
        [0.8,	1.5],
        [0.9,	1.6],
        [1.0,	1.6]
    ]
)

data2 = np.array(
    [
        [0.4,	0.5],
        [0.5,	0.6],
        [0.6,	0.7],
        [0.7,	0.8],
        [0.8,	0.8],
        [0.9,	0.8],
        [1.0,	0.8]
    ]
)

plt.plot(data[:,0], data[:,1], marker='o', label='Sundaram et al 2013', ls='--')
plt.plot(data2[:,0], data2[:,1], marker='o', label='Du et al 2003', ls='--')

plt.grid(which='major', color='green')
plt.minorticks_on()
plt.grid(which='minor', color='green', ls='--')

plt.title('Minimum Ignition Length for Pellet')

plt.xlabel('Particle Volume Fraction')
plt.ylabel('Ignition Length (mm)')

plt.legend()
plt.show()