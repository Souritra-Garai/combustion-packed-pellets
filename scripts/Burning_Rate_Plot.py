import numpy as np
import matplotlib.pyplot as plt

data = [
    [1.0,	40,	21],
    [0.9,	39,	15],
    [0.8,	38,	18],
    [0.7,	36,	13],
    [0.6,	34,	11],
    [0.5,	29,	10],
    [0.4,	26,	55]
]

data_2 = [    
    [0.4,	35,	30],
    [0.5,	49,	28],
    [0.6,	55,	27],
    [0.7,	56,	16],
    [0.8,	58,	20],
    [0.9,	60,	20],
    [1.0,	63,	27]
]

data = np.array(data)
data_2 = np.array(data_2)

plt.errorbar(data[:,0], data[:,1], marker='o', label='Sundaram et al 2013', ls='--', yerr=data[:, 2], uplims=True, lolims=True)
plt.errorbar(data_2[:,0], data_2[:,1], marker='o', label='Du et al 2003', ls='--', yerr=data_2[:, 2], uplims=True, lolims=True)

plt.grid(which='major', color='green')
plt.minorticks_on()
plt.grid(which='minor', color='green', ls='--')

plt.title('Average Velocity of Combustion Front in Pellet')

plt.xlabel('Particle Volume Fraction')
plt.ylabel('Combustion Front Velocity (mm/s)')

plt.legend()
plt.show()