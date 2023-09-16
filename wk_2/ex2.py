import matplotlib
import math
import matplotlib.pyplot as plt
import numpy as np

A = 10
B = 5

NuX = 5
NuY = 17
Alpha = math.pi / 12

def x(t, nu):
    return A * math.sin(nu * t)

def y(t, nu, Alpha):
    return B * math.sin(nu * t + Alpha)

T = np.arange(0, 2 * math.pi, 0.01)

fig, axs = plt.subplots(1, 2)
fig.suptitle('Фигура Лиссажу')

X = [x(t, NuX) for t in T]
Y = [y(t, NuY, Alpha) for t in T]

axs[0].plot(X, Y, 'y')
axs.flat[0].set(xlabel='Alpha = pi / 12', ylabel='5 : 17')
axs[0].grid()

axs[1].hist(X, 50, density=True, facecolor='g', alpha=0.75)
axs.flat[1].set(xlabel='Распределение узловых точек')
axs[1].grid()

plt.show()