import os
import sys

import numpy as np
import matplotlib.pyplot as plt

all_folders = [os.path.join(os.path.dirname(sys.path[0]), 'solutions/' + d) for d in os.listdir(os.path.join(os.path.dirname(sys.path[0]), 'solutions'))]

# print(all_folder)

folder = max(all_folders, key=os.path.getmtime)

print('Visualizing solution from\n' + folder)

T = np.genfromtxt(os.path.join(folder, 'Temperature.csv'), delimiter=',')[:, :-1]
X = np.genfromtxt(os.path.join(folder, 'Conversion.csv'), delimiter=',')[:, :-1]

N, M = T.shape

Dt = 1E-5
Dx = 6.35E-3 / (M-1)
L = 6.35E-3

with open(os.path.join(folder, 'Combustion_Config.txt')) as file :

    for line in file :

        mylist = list(line.split('\t'))

        if mylist[0] == 'Delta t :' :

            Dt = float(mylist[1])

            # print(mylist, Dt)

        if mylist[0] == 'Delta x :' :

            Dx = float(mylist[1])

            # print(mylist, Dt)

        if mylist[0] == 'Length :' :

            L = float(mylist[1])

            # print(mylist, Dt)        

combustion_front = []

for eta in X :

    for i in range(M) :

        if np.isclose(eta[i], 0.000001) :

            break

    combustion_front.append(i*Dx)

combustion_front = np.array(combustion_front)
combustion_velocity = (combustion_front[1:] - combustion_front[:-1]) * 1E3 / Dt

index = np.greater(combustion_velocity, 5, where=np.less(combustion_velocity, 200))

print(np.mean(combustion_velocity[index]))
print(np.std(combustion_velocity[index]))

plt.hist(combustion_velocity, bins=400, range=[5, 200])

plt.grid(which='major', color='green')
plt.minorticks_on()
plt.grid(which='minor', color='green', ls='--')

plt.show()