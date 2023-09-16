import matplotlib
import math
import matplotlib.pyplot as plt
import numpy as np

print("Hello, world!")

A = 10
B = 10

Alpha1 = 0
Alpha2 = math.pi / 4
Alpha3 = math.pi / 2
Alpha4 = 3 * math.pi / 4
Alpha5 = math.pi

Alpha = np.array([Alpha1, Alpha2, Alpha3, Alpha4, Alpha5])
Nus = np.array([[1, 1], [1, 2], [2, 1], [1, 3], [2, 3]])

def x(t, nu):
    return A * math.sin(nu * t)

def y(t, nu, Alpha):
    return B * math.sin(nu * t + Alpha)

T = np.arange(0, 2 * math.pi, 0.01)

fig, axs = plt.subplots(5, 5)
fig.suptitle('Фигуры Лиссажу')

for i in range(5):
    for j in range(5):
        X = [x(t, Nus[i][0]) for t in T]
        Y = [y(t, Nus[i][1], Alpha[j]) for t in T]
        axs[i, j].plot(X, Y, 'y')

axs.flat[0].set(ylabel='1 : 1')
axs.flat[5].set(ylabel='1 : 2')
axs.flat[10].set(ylabel='2 : 1')
axs.flat[15].set(ylabel='1 : 3')
axs.flat[20].set(xlabel='Alpha = 0', ylabel='2 : 3')

axs.flat[21].set(xlabel='Alpha = pi / 4')
axs.flat[22].set(xlabel='Alpha = pi / 2')
axs.flat[23].set(xlabel='Alpha = 3pi / 4')
axs.flat[24].set(xlabel='Alpha = pi')

for ax in axs.flat:
    ax.label_outer()

plt.show()