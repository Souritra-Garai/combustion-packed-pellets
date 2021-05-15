import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

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

fig = plt.figure(figsize=(16,12))
ax = fig.subplots()

t = np.linspace(0, (N-1)*Dt, N)
x = np.linspace(0, L, M)
T_plot, = ax.plot(x, T[0], color='red')
# X_plot, = ax.plot(x, solution.y[100:-1,0]*1000, color='black')

ax.set_xlabel('x (m)')
ax.set_ylabel('Temperature (K)')

ax.set_ylim([0, 3000])

ax.set_title('Temperature Evolution in the Combustion Pellet')

title = ax.text(
    0.5,
    0.9,
    "",
    bbox={'facecolor':'w', 'alpha':1.0, 'pad':5},
    transform=ax.transAxes,
    ha="center"
)

ax.grid(which="major", color='green')
ax.minorticks_on()
ax.grid(which="minor", color='green', ls='--')

def animate(i):

    title.set_text(u"Time t = {:.5f} seconds".format(t[i]))

    T_plot.set_ydata(T[i])  # update the data.
    # X_plot.set_ydata(solution.y[100:-1,i]*1000)

    return T_plot, title, # X_plot

ani = animation.FuncAnimation(
    fig,
    animate,
    interval=50,
    blit=True,
    frames=N
)

# To save the animation, use e.g.
writer = animation.FFMpegWriter(
    fps=100,
    bitrate=1800
)

ani.save("movie.mp4")

# plt.show()